maternal_dynamic <- function(t,y,parms, time.step='month'){
  
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  um=parms$um
  
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  v=parms$VacPro[t]*parms$cover # indicator function for when the vaccine will be administered * vaccine coverage
  RR=parms$RR[t] # relative risk of infection after vaccination; equal to 1 before vaccine introduction
  waningV=parms$waningV
  rrM=parms$rrM[t]
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #??Is this right?? Yes, gamma3 stands for rate of recovery from subsequent infection
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M'] # newborns without vaccination of mothers 
  Mv <-  States[,'Mv'] # newborns with vaccination of mothers 
  Sv <-  States[,'Sv']  # newborns with vaccination of mothers has lower risk of infections when maternal immunity wanes
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  V1 <-  States[,'V1'] # vaccination of pregnant women
  
  N.ages <- length(M)
  
  ###Check the standardization of beta and overall structure of lambda here
  ####################
  seasonal.txn <- (1+parms$b1*cos(2*pi*(t-parms$phi*period)/period))# seasonality waves
  
  # baseline.txn.rate is the probability of transmission given contact per capita
  # (parms$dur.days1/30.44) is the duration of infectiousness of primary infection
  # q depends on transmission type (whether depends on population density or not)
  # density (q=0) vs frequency-dependent (q=1) transmission
  # c2 is the contact matrix 
  # beta is transimissibility per unit time
  transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step)
  beta=transmission_unittime*parms$c2
  
  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)
  infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4
  
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ##########?????????????????##########################  
  
  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  # newborns without vaccination of mothers is still protected from infections because of maternal immunity
  dy[,'M'] <- (1-v)*period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  # newborns with vaccination of mothers is also protected from infections because of maternal immunity
  dy[,'Mv'] <- v*period.birth.rate*sum(States) - 
    (omega+(mu+um))*Mv +
    Aging.Prop*c(0,Mv[1:(N.ages-1)])
  
  # after the maternal immunity wanes, newborns with vaccination of mothers has lower risk of infections 
  dy[,'Sv'] <- omega*Mv - 
    waningV*Sv-
    RR*lambda*Sv -
    (mu + um)*Sv + 
    Aging.Prop*c(0,Sv[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M +
    waningV*Sv-
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0+ RR*lambda*Sv - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    parms$sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- parms$sigma1*lambda*S1 - 
    gamma2*I2-
    (mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 +
    waningV*V1-
    parms$sigma3*lambda*S3 -
    c(rep(0,10),v*period.birth.rate[1]*sum(States),rep(0,2)) -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  # pregnant women (age 20-39) will enter the V1 compartment; the number is equal to newborn (approximately because most of births are single live birth)
  # I did not add any reduced susceptibility or infectiousness at this point and therefore V1 is no different than S3
  dy[,'V1'] <- c(rep(0,10),v*period.birth.rate[1]*sum(States),rep(0,2)) -
    rrM*parms$sigma3*lambda*V1 -
    (mu + um)*V1 -
    waningV*V1+ 
    Aging.Prop*c(0,V1[1:(N.ages-1)]) 
  
  dy[,'I4'] <- parms$sigma3*lambda*S3+rrM*parms$sigma3*lambda*V1 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}