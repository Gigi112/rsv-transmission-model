MAbs_campaign <-  function(t,y,parms, time.step='month'){
  
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
  
  omega = 1/(parms$DurationMatImmunityDays/length.step) # waning rate of maternal immunity
  omega_v=1/(parms$DurationNirImmunityDays/length.step)  # waning rate of monoclonal antibodies
  Vsc=parms$VacPro_campaign[t] #VacPro_campaign is an indicator,
  # all infants under under 6 months will receive a dose in November every year
  Vsr=parms$VacPro_routine[t] #VacPro_routine is an indicator 
  # all infants born between November and March will receive a dose at birth
  Csc=parms$cover # percent of coverage; between 0-1; should be informed by Pediatric vaccine coverage
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #??Is this right?? Yes, gamma3 stands for rate of recovery from subsequent infection
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M']
  Pm <-  States[,'Pm'] 
  P0 <-  States[,'P0'] 
  P1 <-  States[,'P1'] 
  P2 <-  States[,'P2'] 
  P3 <-  States[,'P3'] 
  # infants receive monoclonal antibodies in November if they are under
  # 1 year old
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  N.ages <- length(M)
  
  ##force of infection##################
  
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
  
  # infants under 6 months of age will receive a one time immunization in November each year
  # v is an indicator function for vaccine administrations; will be useful when consider year-round strategy vs seasonal strategy
  # cover indicates the percentage coverage of vaccination; comes from pediatric vaccination coverage data
  dy[,'M'] <- (1-Vsr*Csc)*period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) -
    c(rep(Vsc*Csc,3),rep(0,10))*M
  # a proportion of infants that are 2-3 month old will leave M and S0 and enter P
  # because they receive monoclonal antibodies at that time; correspond to the clinical study
  # https://www.nejm.org/doi/10.1056/NEJMoa1913556?url_ver=Z39.88-2003
  
  # infants who receive monoclonal antibodies at birth will be protected by MAbs
  dy[,'Pm'] <- Vsr*Csc*period.birth.rate*sum(States) - 
    (omega_v+(mu+um))*Pm + # omega_v should be smaller than omega because MAbs has a longer protection period
    Aging.Prop*c(0,Pm[1:(N.ages-1)]) 
  
  # After either of the immunity wanes, infants become susceptible to infections
  dy[,'S0'] <- omega*M+omega_v*Pm+omega_v*P0 -
    lambda*S0 -
    (mu + um)*S0 - 
    c(rep(Vsc*Csc,3),rep(0,10))*S0 +
    Aging.Prop*c(0,S0[1:(N.ages-1)])
  
  
  dy[,'P0'] <-  c(rep(Vsc*Csc,3),rep(0,10))*M +
    c(rep(Vsc*Csc,3),rep(0,10))*S0 - 
    (omega_v+(mu+um))*P0 +
    Aging.Prop*c(0,P0[1:(N.ages-1)]) 
  
  # the rest is the same as the original RSV transmission model
  dy[,'I1'] <-   lambda*S0 - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    parms$sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) +
    omega_v*P1 -
    c(rep(Vsc*Csc,3),rep(0,10))*S1
  
  dy[,'P1'] <-  c(rep(Vsc*Csc,3),rep(0,10))*S1 - 
    (omega_v+(mu+um))*P1 +
    Aging.Prop*c(0,P1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- parms$sigma1*lambda*S1 - 
    gamma2*I2-(mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) +
    omega_v*P2 -
    c(rep(Vsc*Csc,3),rep(0,10))*S2
  
  dy[,'P2'] <-  c(rep(Vsc*Csc,3),rep(0,10))*S2 - 
    (omega_v+(mu+um))*P2 +
    Aging.Prop*c(0,P2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 -
    parms$sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) +
    omega_v*P3 -
    c(rep(Vsc*Csc,3),rep(0,10))*S3
  
  dy[,'P3'] <- c(rep(Vsc*Csc,3),rep(0,10))*S3 - 
    (omega_v+(mu+um))*P3 +
    Aging.Prop*c(0,P3[1:(N.ages-1)]) 
  
  dy[,'I4'] <- parms$sigma3*lambda*S3 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}