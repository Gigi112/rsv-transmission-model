functions {
  real[] msis (real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
               
// state, the volumes in each compartment, y;
// theta, variables used to compute f, which depend on the model parameters;
// x_r, real variables used to evaluate f, which only depend on fixed data;
// x_i, integer values used to evaluate f, which only depend on fixed data.

      real M = y[1];
      real S0 = y[2];
      real I1 = y[3];
      real S1 = y[4];
      real I2 = y[5];
      real S2 = y[6];
      real I3 = y[7];
      real S3 = y[8];
      real I4 = y[9];
      real N = x_i[1];
      
      real birthrate = x_r[1];
      real um = x_r[2];
      real rho1 = x_r[3];
      real rho2 = x_r[4];
      real gamma1 = x_r[5];
      real gamma2 = x_r[6];
      real gamma3 = x_r[7];
      real sigma1 = x_r[8];
      real sigma2 = x_r[9];
      real sigma3 = x_r[10];
      
      real beta = theta[1];
      real omega = theta[2];
      
      real lambda = beta * (I1 + rho1*I2 + rho2*I3 + rho2*I4);
      
      real dM_dt =  birthrate*N - (omega+um)*M;
      real dS0_dt =  omega*M - lambda*S0 - um*S0;
      real dI1_dt =  lambda*S0 - gamma1*I1 - um*I1;
      real dS1_dt =  gamma1*I1 - sigma1*lambda*S1 - um*S1;
      real dI2_dt =  sigma1*lambda*S1 - gamma2*I2 - um*I2;
      real dS2_dt =  gamma2*I2 - sigma2*lambda*S2 - um*S2;
      real dI3_dt =  sigma2*lambda*S2 - gamma3*I3 - um*I3;
      real dS3_dt =  gamma3*I3 + gamma3*I4 - sigma3*lambda*S3 - um*S3;
      real dI4_dt =  sigma3*lambda*S3 - gamma3*I4 - um*I4;
      
      return {dM_dt, dS0_dt, dI1_dt,dS1_dt, dI2_dt, dS2_dt, dI3_dt,dS3_dt, dI4_dt};
  }
}

// Fixed data is declared in the data block:
data {
  int<lower=1> n_days;
  real y0[9];
  real t0;
  real ts[n_days];
  int N;
  real birthrate;
  real um;
  real rho1;
  real rho2;
  real gamma1;
  real gamma2;
  real gamma3;
  real sigma1;
  real sigma2;
  real sigma3;
  int hosp_cases[n_days];
  real hosp1;
  real hosp2;
  real hosp3;
}


transformed data {
  real x_r[10];
  int x_i[1] = { N };
  
  x_r[1]= birthrate;
  x_r[2]=um;
  x_r[3]=rho1;
  x_r[4]=rho2;
  x_r[5]=gamma1;
  x_r[6]=gamma2;
  x_r[7]=gamma3;
  x_r[8]=sigma1;
  x_r[9]=sigma2;
  x_r[10]=sigma3;
}

parameters {
  real<lower=0> beta;
  real<lower=0> omega;
}

transformed parameters{
  real y[n_days, 9]; // y is the states matrix that has row length=n_days and columns=3.
  real lambda[n_days];
  
    // outcomes
  vector[n_days] output_hosp; // overall case incidence by day
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = omega;

    y = integrate_ode_rk45(msis, y0, t0, ts, theta, x_r, x_i);
  }
  
  for(i in 1:n_days){
    lambda[i] = beta * (y[i,3] + rho1*y[i,5] + rho2*y[i,7] + rho2*y[i,9]);
    output_hosp[i] = hosp1*y[i,2]*lambda[i]+hosp2*sigma1*y[i,4]*lambda[i]+hosp3*sigma2*y[i,6]*lambda[i]+hosp3*sigma3*y[i,8]*lambda[i];
  }
}

model {
  //priors
  beta ~ normal(2, 1); //truncated at 0
  omega ~ normal(0.4, 0.5); //truncated at 0
  
  for(i in 1:n_days) {
  //sampling distribution
  target += poisson_lpmf(hosp_cases[i] | output_hosp[i]);
  }
}
