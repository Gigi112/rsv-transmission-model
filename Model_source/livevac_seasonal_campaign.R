livevac_seasonal_campaign<- function(t,y,parms, time.step='month'){
  
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
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #gamma3 stands for rate of recovery from subsequent infection
  gamma_v= 1/(parms$dur.days_v/length.step)
  
  
  Cr = parms$Cr  # vaccine coverage will be informed by pediatric vaccine coverage
  
  # 
  Vrt = parms$VacPro_routine[t]#VacPro_routine is an indicator 
  #  infants will receive a dose at age 2-3 month between November and March
  Vcp = parms$VacPro_campaign[t] #VacPro_campaign is an indicator,
  # all infants under 9 months will receive a dose in November every year
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  Iv <-  States[,'Iv']
  Sv <-  States[,'Sv']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  
  N.ages <- length(M)
  
  ##seasonal lambda parameter##################
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
  lambda_v <- Iv%*%beta_a_i*parms$RR  # lambda_v= transmissibility of vaccine virus * number of vaccinated individuals
  # that shed vaccine virus currently
  ##########?????????????????##########################  
  
  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    c(0,M[1:(N.ages-1)])*Aging.Prop -
    c(0,M[1:(N.ages-1)])*c(0,Cr*Vrt,rep(0,N.ages-2))*Aging.Prop -
    c(0,rep(Vcp*Cr,4),rep(0,8))*M
  
  dy[,'S0'] <- omega*M -
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) -
    c(0,Cr*Vrt,rep(0,N.ages-2))*Aging.Prop*c(0,S0[1:(N.ages-1)]) -
    lambda_v*S0 -
    c(0,rep(Vcp*Cr,4),rep(0,8))*S0
  
  dy[,'I1'] <-   lambda*S0 - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  # Iv stands for just got vaccinated and can shed vaccine virus to others
  # a proportion of infants that are 2-3 month old will leave M and S0 and enter Iv
  # because they receive live attenuated vaccine at that time; correspond to the expected vaccination schedule
  dy[,'Iv'] <-  c(0,M[1:(N.ages-1)])*c(0,Cr*Vrt,rep(0,N.ages-2))*Aging.Prop +
    c(0,Cr*Vrt,rep(0,N.ages-2))*Aging.Prop*c(0,S0[1:(N.ages-1)])  +
    lambda_v*S0 +
    c(0,rep(Vcp*Cr,4),rep(0,8))*M +
    c(0,rep(Vcp*Cr,4),rep(0,8))*S0 -
      ( mu + um)*Iv + 
      Aging.Prop*c(0,Iv[1:(N.ages-1)]) - 
      gamma_v*Iv
  
  # when the infants stop shedding vaccine virus, they enter the Sv compartment
  # we assumed their susceptibility to infections is now equal to after a natural infection (S1) because of vaccination 
  dy[,'Sv'] <-   gamma_v*Iv - 
    parms$sigma1*lambda*Sv -
    (mu+um)*Sv + 
    Aging.Prop*c(0,Sv[1:(N.ages-1)])
  
  # the rest is very similar to the original transmission model
  dy[,'S1'] <- gamma1*I1 - 
    parms$sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)])
  
  dy[,'I2'] <- parms$sigma1*lambda*S1 +
    parms$sigma1*lambda*Sv - 
    gamma2*I2-(mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    parms$sigma2*lambda*S2 -
    (mu+um)*S2+ 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 -
    parms$sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)])
  
  dy[,'I4'] <- parms$sigma3*lambda*S3 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  
  
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}
