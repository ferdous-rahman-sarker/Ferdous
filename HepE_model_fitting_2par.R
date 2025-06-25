library(deSolve)

cases <- read.csv('hepEdata_begin.csv')

hepENLL<- function(par)
{

  # define the number of weeks to run the model
  times <- seq(0, length(cases$week), by = 1)
  week_start <- 38
  week_stop <- 78
  #MODEL PARAMETERS: ALL
  parameters <- c(mui=(1/(50*52)),    # birth
                  muo=(1/(50*52)),    # death
                  R0=exp(par[1]),     # basic reproduction number
                  omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                  gamma=(1/2),        # rate of movement from latent to infectious stage = 1/(average latent period)
                  tau=(1/4),          # rate of recovery = 1/(average duration of infection)
                  report=1/7,         # proportion of all infections that are reported
                  amp=0,              # relative amplitude of seasonal forcing
                  phi=0              # week of peak in seasonal forcing
  )
  
  # MODEL INITIAL CONDITIONS
  initP<-exp(par[2])
  
  initE<-1
  initI<-0
  initR<-0
  initS<-initP-initE-initI-initR
  
  state <- c(S = initS, E=initE, I = initI,R = initR)
  
  # set up a function to solve and fit the model
  HepE<-function(t, state, parameters) 
  {
    with(as.list(c(state, parameters)),
{
  
  # define variables
  P <- (S+E+I+R)
  seas<-1+amp*cos(2*pi*(t-phi)/52)
  beta<-(R0*(muo+tau)*(gamma+muo))/gamma
  lam <- beta*I/P
  
  # rate of change
  dS <- mui*P-muo*S-lam*S+omega*R
  dE <- -muo*E+lam*S-gamma*E
  dI <- -muo*I+gamma*E-tau*I
  dR <- -muo*R+tau*I-omega*R
  
  # return the rate of change
  list(c(dS, dE, dI, dR))
}
    ) 
    
  }  
  
  out <- ode(y = state, times = times, func = HepE, parms = parameters)
  pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]

  
  beta  <- (parameters["R0"]*(parameters["muo"]+parameters["tau"])*(parameters["gamma"]+parameters["muo"]))/parameters["gamma"]
  lam <- beta*out[,"I"]/pop
  inc <- parameters["report"]*lam*out[,"S"]
  
  plot(cases$week,cases$new.case,pch=19,col='red',xlim=c(week_start,week_stop),ylim=c(0,120),main = "Hepatitis E Outbreak: model versus data",xlab = "Time",ylab="new reported cases")
  lines(out[,"time"]+week_start,inc,lwd=3)
  
  Sys.sleep(.1)
  
  L=matrix(0,nrow=length(cases[,1]),ncol=1)  
  L<-dpois(cases$new.case,inc,log=TRUE)
  
  LL<-sum((1-(L==-Inf))*L,na.rm=TRUE)
  
  NLL<--LL
  
  NLL 
  
}

fit<-optim(c(log(4),log(209100)), hepENLL, hessian=TRUE,method=("Nelder-Mead"),control=list(trace=TRUE))

fisher_info <- solve(fit$hessian)
prop_sigma <- sqrt(diag(fisher_info))
upper <- exp(fit$par+1.96*prop_sigma)
lower <- exp(fit$par-1.96*prop_sigma)
estimate<-exp(fit$par)
estimate
