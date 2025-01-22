library(jagsUI)

######################
# simulated count data
######################

n.sites = 100 # scalar, number of habitat patches
n.years = 10 # scalar, number of years
nreps = 3 # scalar, number of surveys within a primary sampling period (needed to estimate observation error)



#################################

# storage arrays
N = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the true counts of adults
S = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the surviving adults from previous year
R = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits

H = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of fished animals
IH = array(dim = c(n.sites, n.years))
r = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits


# parameters

lambda = 100          # mean number of adults in year 1 
psi = 0.9             # adult survival
gamma = 40            # recruitment
F = 0.5               # total fishing pressure
p = 0.9               # detection
alpha = 0.03          # reporting ~ N slope



# create simulated data

# year 1
N[,1] = rpois(n.sites, lambda) # if pop sizes across sites are overdispersed, can't use poisson

for(t in 1:(n.years-1)){
  
  R[,t+1] = rpois(n.sites, gamma)
  

  for(s in 1:n.sites){
    
    psi.tmp = rnorm(1, psi, 0.01)
    F.tmp = rnorm(1, F, 0.05)
    A.tmp = rnorm(1, alpha, 0.005)
    
    r[s,t] = 1/(1+exp(alpha*-N[s,t]))

    psi.tmp = rnorm(1, psi, 0.01)
    F.tmp = rnorm(1, F, 0.01)
    ps = c(psi.tmp-F.tmp, F.tmp*r[s,t], F.tmp*(1-r[s,t]), (1-psi.tmp))
    ps = ps/sum(ps)

  m = rmultinom(1, N[s,t], prob = ps)

      
    S[s,t+1] = m[1] 
    H[s,t+1] = m[2]
    IH[s,t+1] = m[3]
    
    N[s,t+1] = (S[s,t+1] + R[s,t+1])
    
  }
}

# true declines
IHprop.TRUE =  N[,1] - N[,10]
IHprop.TRUE = ifelse(IHprop.TRUE < 0, 0, 1)
sum(IHprop.TRUE)

r.mean = mean(r, na.rm = T)
y = array(dim = c(n.sites, n.years, nreps))

ps
mean(N)
mean(S, na.rm = T) / mean(N)
median(H, na.rm = T) / mean(N)
median(IH, na.rm = T) / mean(N)


# simulate observation error of counts (but to accurately estimate, you need repeated counts within a year)
for(j in 1:nreps){
  y[,,j] = rbinom(n.sites*n.years, N, p)     # observed counts of adults
}

plot(N, y[,,1]) # observed vs true counts of adults
abline(0,1)

# observed catch
C = round(H)


op <- par(mfrow = c(2, 2), mar = c(5, 5, 4, 3), cex.lab = 1.5, 
          cex.axis = 1.5)
on.exit(par(op))
matplot(t(N), type = "l", main = paste("Population trajectories"), 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "N")
matplot(t(S), type = "l", main = "Number of apparent survivors", 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "S")
hist(N[, 1], main = "Distribution of N in first year", 
     breaks = 50, col = "grey")
hist(N[, n.years], main = "Distribution of N in last year", 
     breaks = 50, col = "grey")


op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 3), cex.lab = 1.5, 
          cex.axis = 1.5)
on.exit(par(op))
matplot(t(N), type = "l", main = paste("Population trajectories"), 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "N", ylim = c(0, 100))
matplot(t(C), type = "l", main = paste("Reported Catch"), 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "C", ylim = c(0, 100))

############
# Model Runs
############

#Bundle data
jags.data <- list(y = y,                                        # counts of adults and juveniles
                  C=C, gamma = gamma, lambda = lambda, psi = psi,                                        # catch data
                  n.years=ncol(y), n.sites = nrow(y), nreps = nreps) 


#Initial values
Rst <- apply(y, c(1,2), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA

Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N[,1]

inits <- function(){list( 
                         alpha = 40, F = 0.6,  
                         p = p, R = Rst+1, N = Nst+2)}


#parameters monitored
parameters <- c("F", "p", "N")

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 100 ; nb <- 20000 ; nc <- 3


#################
# baseline
###################

#call JAGS from R
(out= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "commercial_baseline.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out)
save(out, file = "commercial_baseline.RData")
load("commercial_baseline.RData")


# posteriors
lest = out$sims.list$lambda # true value 20
pest = out$sims.list$p # true value 0.9
fest = out$sims.list$F # true value 0.5
psiest = out$sims.list$psi # true value 0.8


png("commercial_baseline.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(2,2))

hist(pest, breaks = 10, 
     main = "Detection", xlim = c(0.5, 1))
abline(v = 0.9, col = "red", lwd = 2)

hist(psiest, breaks = 10, 
     main = "Survival", xlim = c(0, 1))
abline(v = 0.8, col = "red", lwd = 2)

hist(fest, breaks = 10, 
     main = "Total Fishing", xlim = c(0, 1))
abline(v = 0.5, col = "red", lwd = 2)

hist(lest, breaks = 10, 
     main = "Initial Pop Size", xlim = c(10, 30))
abline(v = 20, col = "red", lwd = 2)

dev.off()


#######################################
# model accounts for unreported harvest
########################################

#parameters monitored
parameters <- c("F", "alpha", "p", "N", "r")

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 100 ; nb <- 20000 ; nc <- 3

#call JAGS from R
(out2 = jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "commercial_with_unreported.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out2)
save(out2, file = "commercial_with_unreported.RData")
load("commercial_with_unreported.RData")

# posteriors
#lest = out2$sims.list$lambda # true value 20
pest = out2$sims.list$p # true value 0.9
fest = out2$sims.list$F # true value 0.5
#psiest = out2$sims.list$psi # true value 0.8




hist(pest, breaks = 10, 
     main = "Detection", xlim = c(0.5, 1))
abline(v = p, col = "red", lwd = 2)



# posteriors
aest = out2$sims.list$alpha # true value 0.05
rest = out2$sims.list$r # true value 0.9

png("commercial_with_unreported.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(1,2))

hist(aest, breaks = 10, 
     main = "alpha", xlim = c(0, 50))
abline(v = 1/alpha, col = "red", lwd = 2)

hist(fest, breaks = 10, 
     main = "Total Fishing", xlim = c(0, 1))
abline(v = F, col = "red", lwd = 2)

dev.off()


plot(aest, fest)

hist(rest, breaks = 10, 
     main = "reporting", xlim = c(0, 1))
abline(v = r.mean, col = "red", lwd = 2)


#######
# r ~ N
library(scales)

png("commercial_reporting_predictions.png", width = 7, height = 5, units = "in", res = 600)

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

x = 10:100

plot(x, rep(-10, length(x)), type = "p", ylab = "", xlab = " ", cex = 1.5,
     ylim = c(0.5, 1), xlim = c(0, 100), lwd = 2, pch = 5, axes = FALSE, main = " ")

axis(1)
mtext("Abundance", side = 1, line = 3, cex = 1.5, font = 2)
axis(2)
par(las = 0)
mtext("Reporting Rate", side = 2, line = 4, cex = 1.5, font = 2)

lines(x, 1/(1+exp(alpha*-(x))), lwd =2, col = alpha("black", 0.8))
lines(x, 1/(1+exp(-(x)/median(aest))), lwd = 2, col = alpha("dark green", 0.8))
lines(x, 1/(1+exp(-(x)/max(aest))), lwd = 2, col = alpha("dark green", 0.5), lty = "dashed")
lines(x, 1/(1+exp(-(x)/min(aest))), lwd = 2, col = alpha("dark green", 0.5), lty = "dashed")

legend("topleft", legend=c("True Relationship","Predicted Relationship"), col=c("black", "dark green"), pt.cex=2, pch=15)

dev.off()


rs = out2$sims.list$r # true value 0.9
ns = out2$sims.list$N # true value 0.9

df = data.frame(ns = tidyr::gather(as.data.frame(ns[,,1:9]))[,2], rs = tidyr::gather(as.data.frame(rest))[,2])
## Use densCols() output to get density at each point
x <- densCols(df$ns,df$rs, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#FF3100", "#FF9400", "#FCFF00", 
                            "#45FE4F", "#00FEFF", "#000099"))(6)
df$col <- ifelse(df$dens >= 250, cols[1], ifelse(df$dens >= 200, cols[2], ifelse(df$dens >= 150, cols[3], ifelse(df$dens >= 100, cols[4], ifelse(df$dens >= 50, cols[5], cols[6])))))

## Plot it, reordering rows so that densest points are plotted on top
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(rs~ns, data=df[order(df$dens),], pch=20, col=col, cex=1, ylim = c(0.5, 1), 
     xlab = "Abundance", ylab = "Reporting Rate")
lines(40:100, 1/(1+exp(alpha*-(40:100))), lwd = 2)

# estimate bandwidths
h <- c(KernSmooth::dpik(df$ns), KernSmooth::dpik(df$rs))

# obtain density estimatte
f1 <- KernSmooth::bkde2D(df[, 1:2], bandwidth = h, gridsize = c(10000, 10000))

# setup vector of density levels to obtain contours for
contour_levels <- pretty(f1$fhat, 9)

# setup color palette
crp <- colorRampPalette(rev(c("#FF3100", "#FF9400", "#FCFF00", 
                              "#45FE4F", "#00FEFF", "#000099")))

# density based colors 
df$col <- densCols(df$ns, df$rs, bandwidth = h, colramp = crp)

# plot
plot(rs ~ ns, data = df, pch = 20, col = col, cex = 0.5)
contour(f1$x1, f1$x2, f1$fhat, levels = contour_levels, col = crp(9),
        lwd = 2, add = T) 
#############
# comparison
############

#############
# comparison
############

load("commercial_baseline.RData")
load("commercial_with_unreported.RData")

# posteriors
fest = out$sims.list$F # true value 0.5
rest = 1

# posteriors
fest2 = out2$sims.list$F # true value 0.5
rest2 = out2$sims.list$r # true value 0.8

png("commercial_comparison.png", width = 12, height = 8, units = "in", res = 600)

par(mfrow = c(2,3))
op <- par(cex.main = 1.5, mar = c(3, 6, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.2, bty = "n", las = 1)


hist(r, breaks = 20, col=rgb(1,1,1,0.5), 
     main = "Reporting Rate \n True", xlim = c(0, 1), xlab = "", ylab = "")

hist(c(rep(0.99, 20500),rep(1, 20500)), breaks = 1, col=rgb(1,0,0,0.5), 
     main = "Reporting Rate \n assuming 100% reporting", xlim = c(0, 1), ylim = c(0,41000), xlab = "", ylab = "")

hist(rest2, breaks = 20, col=rgb(0,0,1,0.5), 
     main = "Reporting Rate \n assuming some unreporting", xlim = c(0, 1), xlab = "", ylab = "")

hist(0, breaks = 10, col=rgb(1,1,1,0.5), 
     main = "Total Fishing \n True", xlim = c(0.4, 0.6), xlab = "")
abline(v = F, col = "red", lwd = 2)

hist(fest, breaks = 10, col=rgb(1,0,0,0.5), 
     main = "Total Fishing \n assuming 100% reporting", xlim = c(0.4, 0.6), xlab = "", ylab = "")

hist(fest2, breaks = 10, col=rgb(0,0,1,0.5), 
     main = "Total Fishing \n assuming some unreporting", xlim = c(0.4, 0.6), xlab = "", ylab = "")

dev.off()





#######################################
# model misspecification runs
########################################

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 100 ; nb <- 20000 ; nc <- 3

#parameters monitored
parameters <- c("F", "p", "N", "r")

#call JAGS from R
(out.com.sub = jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "subsistence_with_unreported.txt", 
             n.iter=ni, n.burnin = nb, n.chains= nc,
             n.thin= nt, n.adapt = na, parallel = TRUE))

save(out.com.sub, file = "commercial_mis_subsistence.RData")


load("commercial_mis_subsistence.RData")

cms.mu = out.com.sub$q50$F - 0.5
cms.upr = out.com.sub$q2.5$F - 0.5
cms.lwr = out.com.sub$q97.5$F - 0.5

