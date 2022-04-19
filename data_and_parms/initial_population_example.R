# https://data.census.gov/cedsci/table?t=Age%20and%20Sex&g=0400000US34&tid=ACSST5Y2010.S0101
agep=c(rep(6.3/5/6,6),6.3/5,6.3/5*3,6.5,6.7+6.9,6.0+6.3+6.3+7.1,7.7+8.1+7.4+6.2,5.1+3.8+2.9+2.6+2.1+1.9)/100

N=7365011*agep

N_ages <- length(agep)
StateNames <- c('M','S0','I1','S1','I2','S2','I3','S3','I4')
States <- array(NA, dim=c(N_ages, length(StateNames) )) #  N age groups xN parameters 
dimnames(States)[[1]] <- agenames
dimnames(States)[[2]] <- StateNames

yinit.matrix <- array(NA, dim=c(N_ages, length(StateNames) ))

dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

yinit.matrix[,c('S1','I2','S2','I3','S3','I4')]  = 0
yinit.matrix[,'M'] = c(N[1:6], rep(0,N_ages-6))
yinit.matrix[,'S0'] = c(rep(0,6),N[7:N_ages]-rep(N_ages-6)) 
yinit.matrix[,'I1'] = c(rep(0,6), rep(1,N_ages-6))  #initializes with 1 infected person per age group 

yinit.vector <- as.vector(yinit.matrix) #Vectorize the ynit matrix