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
      
      return {lambda,dM_dt, dS0_dt, dI1_dt,dS1_dt, dI2_dt, dS2_dt, dI3_dt,dS3_dt, dI4_dt};
  }
}
