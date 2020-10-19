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


set.seed(20)
n=200;

resolution=n #means we compute at every 1/n

g=(1:resolution)/resolution

a0 <- 25*exp(-(g-0.5)^2/0.1)
a1 <- 0.3*(g-1)^2+0.1
b1 <- 0.1*(g)^1.5+0.1

e=rnorm(n);x=e;sigma2=array(0,n)
sigma20=a0[1]/(1-b1[1]);x0=rpois(1, sigma20);
sigma2[1]=a0[1]+a1[1]*x0+b1[1]*sigma20;x[1]=rpois(1, sigma2[1])

for (i in 2:n)
{
  sigma2[i]=a0[i]+a1[i]*x[i-1]+b1[i]*sigma2[i-1]
  x[i]= rpois(1, sigma2[i])
}  

fitS <- fit.tvINGARCHMCMCcombo(as.numeric(c(x0,x)), order1 = 1, order2 = 1, norder = 4, knot = 6)

mufn_p  <- matrix(unlist(fitS$Mfn), length(fitS$Mfn[[1]]))
Afn_p   <- matrix(unlist(fitS$Afn), length(fitS$Afn[[1]]))
Bfn_p   <- matrix(unlist(fitS$Bfn), length(fitS$Bfn[[1]]))

mufn_p  <- matrix(unlist(fitS$Mfn), length(fitS$Mfn[[1]]))
Afn_p   <- matrix(unlist(fitS$Afn), length(fitS$Afn[[1]]))
Bfn_p   <- matrix(unlist(fitS$Bfn), length(fitS$Bfn[[1]]))

mutp <- rowMeans(mufn_p) #Posterior mean of mean function
Atp  <- rowMeans(Afn_p)  #Posterior mean of AR(1)
Btp  <- rowMeans(Bfn_p)  #Posterior mean of CH(1)
