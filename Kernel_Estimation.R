# Import library
library(pracma) 
library(tidyr)
library(ggplot2)

SP500 <- read.csv(file.choose(), header = TRUE)
head(SP500)

#Import data
p <- (SP500$close)
for (i in 1:(length(p)-1)){
  if (p[i]==p[i+1]){
    p[i] <- (p[i-1]+p[i+1])/2
  }
}

inc <- diff(log(p))
inc <- inc^2
x <- log(inc)

#Create characteristic function 
fi<- function(t=real){
  fi<- (1/sqrt(pi))*2^{1i*t}*gammaz(0.5 + t*1i)
}

#Create Kernel function
fiker<-function(t=real, h=real){
  if ((h*t>=-1)&(h*t<1))
  {fiker<-(1-(h*t)^2)^3} else {fiker<-0}
}

#data discretization
n<-length(x)
M <- 2^7
h <- 0.3 
a <- 2*max(abs(x))

delta <- a/M

t <- array(0,c(M,1))
for ( k in 1:M ){
  t[k] <- (k-1-(M/2))*delta
}

weight <- array(0, c(M,1))
for (k in 1:(M-1)){
  for (i in 1:n){
    if (x[i] >= t[k] & x[i] <= t[k+1]) {
      weight[k] <- weight[k] + ((1/n)*(1/delta^2)*(t[k+1]-x[i]))
      weight[k+1] <- weight[k+1] + ((1/n)*(1/delta^2)*(x[i]-t[k]))
    }
  }
}

#check
sum(weight)
1/delta


# compute empirical characteristic function
for (i in 1:length(weight)){
  if (i %% 2 == 0)
  {weight[i] <- weight[i]} else {weight[i] <- (-1)*weight[i]}
}

Y_l <- (fft(weight, inverse=TRUE))/M

# compute zeta_l*
s_l <- array(0,c(M,1))
for (l in 1:length(s_l)){
  s_l[l] <- (2*pi*(l-1-(M/2)))/a
}

zeta_l <- array(0,c(length(s_l),1))
for (i in 1:length(s_l)){
  zeta_l[i] <- ((fiker(s_l[i],h)/fi(s_l[i])))*Y_l[i]
}

# compute zeta
f_x <-Re(fft(zeta_l, inverse = FALSE))
for (i in 1:length(f_x)){
  if (i %% 2 == 0)
  {f_x[i] <- f_x[i]} else {f_x[i] <- (-1)*f_x[i]}
}
f_x.real <-data.frame(f_x)


#  Plot return's volatility by Kernel method
qplot(t,f_x,geom="line",
      colour = I("black"),
      fill = I("gray44"),
      size=I(1),
      ylab = "density",
      xlab = "")


