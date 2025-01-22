#########################################################
# commerical fishing (model assumes no unreported harvest)

cat(file="commercial_baseline.txt", " 
    model {
    #priors

##########
phi_raw[1] <- psi * (psi-F)
phi_raw[2] <- psi * F
phi_raw[3] <- 0
phi_raw[4] ~ dgamma(0.1, 0.1)

for (i in 1:4) {
  phi[i] <- phi_raw[i] / sum(phi_raw)
}

    
    F ~ dnorm(0.6, 1)I(0,1)

    p ~ dnorm(0.9, 10)I(0,1)
    
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  

    # year 1
  N[s,1] ~ dpois(lambda) 

  for(t in 1:(n.years-1)){
    
    S[s,t+1] ~ dbin(phi[1], N[s,t])
    C[s,t+1] ~ dbin(phi[2], N[s,t]) 
    D[s,t+1] <- N[s,t] - (S[s,t+1] + C[s,t+1]) 
    
    R[s,t+1] ~ dpois(gamma)

    N[s,t+1] <- (S[s,t+1] + R[s,t+1])

    } # t
    

     for (t in 1:n.years) {
      for (j in 1:nreps) {

    y[s,t,j] ~ dbin(p, N[s,t])

      } # j
     } # t
    } # s
}  
    ")



########################################################################################
# commerical fishing (model assumes some unreported harvest that changes with abundance)

cat(file="commercial_with_unreported.txt", " 
    model {
    #priors

    p ~ dbeta(20, 2)

    alpha ~ dgamma(1,1)  

    F ~ dbeta(10,10)I(0,psi)

# model for the initial population size: poisson priors

for (s in 1:n.sites){  

   # year 1
   N[s,1] ~ dpois(lambda) 

  for(t in 1:(n.years-1)){
  
#######################################
# fishing pressure changes with abundance

    r[s,t] = 1/(1+exp(-(N[s,t]/alpha)))

    phi[s,t,1] <- psi-F
    phi[s,t,2] <- F*r[s,t]
    phi[s,t,3] <- F*(1-r[s,t])
    phi[s,t,4] <- 1-psi

#######################################
    
    S[s,t+1] ~ dbin(phi[s,t,1], N[s,t])
    C[s,t+1] ~ dbin(phi[s,t,2]/(1 - phi[s,t,1]), N[s,t]-S[s,t+1])
    UC[s,t+1] ~ dbin(phi[s,t,3]/(1 - (phi[s,t,1] + phi[s,t,2])), N[s,t]-(S[s,t+1]+C[s,t+1])) 
    D[s,t+1] <- N[s,t] - (S[s,t+1] + C[s,t+1] + UC[s, t+1]) 

    R[s,t+1] ~ dpois(gamma)

    N[s,t+1] <- (S[s,t+1] + R[s,t+1])

    } # t
    

     for (t in 1:n.years) {
      for (j in 1:nreps) {

    y[s,t,j] ~ dbin(p, N[s,t])

      } # j
     } # t
    } # s
}  
    ")

