set.seed(1)

Ti <- 1000
t <- 0:Ti/Ti
mut0 <- exp(-(t-0.5)^2/0.1)
at10 <- 0.3*(t-1)^2 + 0.1
at20 <- rep(0, Ti+1) #0.5*t
datavect <- rep(0, Ti+1)
datavect[1] <- rpois(1, mut0[1])
datavect[2] <- rpois(1, mut0[2]+at10[2]*datavect[1])
for(i in 3:Ti){
  datavect[i] <- rpois(1, mut0[i]+at10[i]*datavect[i-1]+at20[i]*datavect[i-2])
}

data <- datavect[-1]

fit      <- fit.tvARMCMCunbd(data, order= 1)
Afnp     <- Reduce('+', fit$Afn) / length(fit$Afn)       #Estimated AR coef fns
mufnp    <- Reduce('+', fit$Mfn) / length(fit$Mfn)       #Estimated mean fn
mufnpder <- Reduce('+', fit$Mfnder) / length(fit$Mfnder) #Estimated derivative of mean fn
Afnpder  <- Reduce('+', fit$Afnder) / length(fit$Afnder) #Estimated derivative of AR coef fns