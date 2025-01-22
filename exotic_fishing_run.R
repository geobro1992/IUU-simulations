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

H = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of legally fished animals
IH = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of illegaly fished animals
r = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits
FT = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits
FR = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits

# parameters

lambda = 100        # mean number of adults in year 1 
psi = 0.9           # adult survival
gamma = 40          # recruitment
F = 2               # total fishing pressure
p = 0.9             # detection
aT = -0.03          # total fishing ~ N slope
alpha = -0.05       # reporting ~ N slope


# create simulated data

# year 1
N[,1] = rpois(n.sites, lambda) # if pop sizes across sites are overdispersed, can't use poisson

for(t in 1:(n.years-1)){
  
  R[,t+1] = rpois(n.sites, gamma)
  
  
  for(s in 1:n.sites){
    
    psi.tmp = rnorm(1, psi, 0.01)
    A.tmp = rnorm(1, alpha, 0.005)
    F.tmp = rnorm(1, F, 0.05)
    
    FT[s,t] = F.tmp/(1+exp(aT*-N[s,t]))
    FR[s,t] = F.tmp/(1+exp(alpha*-N[s,t]))
    
    r[s,t] = 1-(FR[s,t]/FT[s,t])
    
    psi.tmp = rnorm(1, psi, 0.01)

    ps = c(psi.tmp-FT[s,t], FT[s,t]*r[s,t], FT[s,t]*(1-r[s,t]), (1-psi.tmp))
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

hist(r)
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

############
# Model Runs
############

#Bundle data
jags.data <- list(y = y,                                        # counts of adults and juveniles
                  C=C, gamma = gamma, lambda = lambda, psi = psi, aT = aT,                                        # catch data
                  n.years=ncol(y), n.sites = nrow(y), nreps = nreps) 


#Initial values
Rst <- apply(y, c(1,2), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA

Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N[,1]

inits <- function(){list( 
  F = 0.5,  
  p = p, R = Rst+1, N = Nst+2)}


#parameters monitored
parameters <- c("F", "p", "N")

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3


#################
# baseline
###################

#call JAGS from R
(out= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "exotic_baseline.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out)
save(out, file = "exotic_baseline.RData")


#######################################
# model accounts for unreported harvest
########################################

inits <- function(){list( 
  alpha = 0.05, aT = 0.03,  
  p = p, R = Rst+1, N = Nst+2)}

#Bundle data
jags.data <- list(y = y,                                        # counts of adults and juveniles
                  C=C, gamma = gamma, lambda = lambda, psi = psi, Fmax = F,
                  n.years=ncol(y), n.sites = nrow(y), nreps = nreps) 

#parameters monitored
parameters <- c("aT", "F", "alpha", "p", "N", "r")

#MCMC settings
na <- 10000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3

#call JAGS from R
(out2 = jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "exotic_with_unreported.txt", 
             n.iter=ni, n.burnin = nb, n.chains= nc,
             n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out2)
save(out2, file = "exotic_with_unreported.RData")
load("exotic_with_unreported.RData")

# posteriors
pest = out2$sims.list$p 
fest = out2$sims.list$F
aest = out2$sims.list$alpha 
aTest = out2$sims.list$aT 
rest = out2$sims.list$r 



hist(pest, breaks = 10, 
     main = "Detection", xlim = c(0.5, 1))
abline(v = p, col = "red", lwd = 2)




par(mfrow = c(1,2))
hist(r, breaks = 10, 
     main = "true reporting", xlim = c(0, 1))

hist(rest, breaks = 10, 
     main = "predicted reporting", xlim = c(0, 1))


#######
# r ~ N


x = 10:100

r.true = 1-(F/(1+exp(alpha*-(x))))/(F/(1+exp(aT*-(x))))

r.pred = 1-(F/(1+exp(-median(aest)*-(x))))/(F/(1+exp(-median(aTest)*-(x))))
r.predU = 1-(F/(1+exp(-max(aest)*-(x))))/(F/(1+exp(-max(aTest)*-(x))))
r.predL = 1-(F/(1+exp(-min(aest)*-(x))))/(F/(1+exp(-min(aTest)*-(x))))


png("exotic_reporting_predictions.png", width = 7, height = 5, units = "in", res = 600)

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)


plot(x, rep(-10, length(x)), type = "p", ylab = "", xlab = " ", cex = 1.5, 
     ylim = c(0, 1), xlim = c(0, 100), lwd = 2, pch = 5, axes = FALSE, main = " ")

axis(1)
mtext("Abundance", side = 1, line = 3, cex = 1.5, font = 2)
axis(2)
par(las = 0)
mtext("Reporting Rate", side = 2, line = 4, cex = 1.5, font = 2)

lines(x, r.true, lwd = 2, col = alpha("black", 0.8))
lines(x, r.pred, lwd = 2, col = alpha("dark green", 0.8))
lines(x, r.predU, lwd = 2, col = alpha("dark green", 0.5), lty = "dashed")
lines(x, r.predL, lwd = 2, col = alpha("dark green", 0.5), lty = "dashed")

legend("topleft", legend=c("True Relationship","Predicted Relationship"), col=c("black", "dark green"), pt.cex=2, pch=15)

dev.off()

#############
# comparison
############

load("exotic_baseline.RData")
load("exotic_with_unreported.RData")

# posteriors
fest = out$sims.list$F # true value 0.5
rest = 1

# posteriors
fest2 = out2$sims.list$F # true value 0.5
rest2 = out2$sims.list$r # true value 0.8

png("exotic_comparison.png", width = 12, height = 8, units = "in", res = 600)

par(mfrow = c(2,3))
op <- par(cex.main = 1.5, mar = c(3, 6, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.2, bty = "n", las = 1)


hist(r, breaks = 20, col=rgb(1,1,1,0.5), 
     main = "Reporting Rate \n True", xlim = c(0, 1), xlab = "", ylab = "")

hist(c(rep(0.99, 250000),rep(1, 250000)), breaks = 1, col=rgb(1,0,0,0.5), 
     main = "Reporting Rate \n assuming 100% reporting", xlim = c(0, 1), ylim = c(0,500000), xlab = "", ylab = "")

hist(rest2, breaks = 20, col=rgb(0,0,1,0.5), 
     main = "Reporting Rate \n assuming some unreporting", xlim = c(0, 1), xlab = "", ylab = "")

hist(FT, breaks = 10, col=rgb(1,1,1,0.5), 
     main = "Total Fishing \n True", xlim = c(0, 1), xlab = "")

hist(fest, breaks = 1, col=rgb(1,0,0,0.5), 
     main = "Total Fishing \n assuming 100% reporting", xlim = c(0, 1), xlab = "", ylab = "")

hist(fest2, breaks = 10, col=rgb(0,0,1,0.5), 
     main = "Total Fishing \n assuming some unreporting", xlim = c(0, 1), xlab = "", ylab = "")
abline(v = F, col = "red", lwd = 2)

dev.off()



#######################################
# model misspecification runs
########################################

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 100 ; nb <- 20000 ; nc <- 3

#parameters monitored
parameters <- c("F", "p", "N", "r")

#call JAGS from R
(out.exo.sub = jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "subsistence_with_unreported.txt", 
                    n.iter=ni, n.burnin = nb, n.chains= nc,
                    n.thin= nt, n.adapt = na, parallel = TRUE))

save(out.exo.sub, file = "exotic_mis_subsistence.RData")


load("exotic_mis_subsistence.RData")

F.true = mean(FT[,-10])

ems.mu = out.exo.sub$q50$F - F.true
ems.upr = out.exo.sub$q2.5$F - F.true
ems.lwr = out.exo.sub$q97.5$F - F.true
