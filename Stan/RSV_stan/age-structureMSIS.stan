functions {
  real[] msis (real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {
               
// state, the volumes in each compartment, y;
// theta, variables used to compute f, which depend on the model parameters;
// x_r, real variables used to evaluate f, which only depend on fixed data;
// x_i, integer values used to evaluate f, which only depend on fixed data.
      int agegroups = x_i[1];
      int q=x_i[2];
      
      real birthrate[agegroups] = x_r[1:agegroups];
      real um = x_r[agegroups+1];
      real rho1 = x_r[agegroups+2];
      real rho2 = x_r[agegroups+3];
      real gamma1 = x_r[agegroups+4];
      real gamma2 = x_r[agegroups+5];
      real gamma3 = x_r[agegroups+6];
      real sigma1 = x_r[agegroups+7];
      real sigma2 = x_r[agegroups+8];
      real sigma3 = x_r[agegroups+9];
      real omega = x_r[agegroups+10];
      real u[agegroups]= x_r[(agegroups+11):(11+2*agegroups)];
      real c2[agegroups*agegroups]=x_r[(12+2*agegroups):(12+2*agegroups+agegroups*agegroups)];
      
      real beta = theta[1];
      real b1 = theta[2];
      real trans = theta[3];
      real M[agegroups] = y[1:agegroups];
      real S0[agegroups] = y[(agegroups+1):(2*agegroups)];
      real I1[agegroups] = y[(2*agegroups+1):(3*agegroups)];
      real S1[agegroups] = y[(3*agegroups+1):(4*agegroups)];
      real I2[agegroups] = y[(4*agegroups+1):(5*agegroups)];
      real S2[agegroups] = y[(5*agegroups+1):(6*agegroups)];
      real I3[agegroups] = y[(6*agegroups+1):(7*agegroups)];
      real S3[agegroups] = y[(7*agegroups+1):(8*agegroups)];
      real I4[agegroups] = y[(8*agegroups+1):(9*agegroups)];
      
      real dydt[9*agegroups];
      real lambda[agegroups];
      real InfectN[agegroups];
       
      real phi = (2*pi()*exp(trans))/(1+exp(trans));
      real season_txn = (1+b1*cos(2*pi()*t/12 - phi));
      
      real beta_contact = beta/(sum(y)^(1-q));
      for (k in 1:agegroups) {
      InfectN[k] = (I1[k]+rho1*I2[k]+rho2*I3[k]+rho2*I4[k]);
      } 
      
      for (a in 1:agegroups) {
      
       lambda[a] = season_txn*sum(to_vector(c2[((a-1)*agegroups+1):(a*agegroups)]).*to_vector(InfectN));
      
       dydt[a] =  log(birthrate[a]+1)/12*sum(y) - (omega+u[a]+um)*M[a];
       dydt[a+agegroups] =  omega*M[a] - lambda[a]*S0[a] - (um+u[a])*S0[a];
       dydt[a+2*agegroups] =  lambda[a]*S0[a] - gamma1*I1[a] - (um+u[a])*I1[a];
       dydt[a+3*agegroups] =  gamma1*I1[a] - sigma1*lambda[a]*S1[a] - (um+u[a])*S1[a];
       dydt[a+4*agegroups] =  sigma1*lambda[a]*S1[a] - gamma2*I2[a] - (um+u[a])*I2[a];
       dydt[a+5*agegroups] =  gamma2*I2[a] - sigma2*lambda[a]*S2[a] - (um+u[a])*S2[a];
       dydt[a+6*agegroups] =  sigma2*lambda[a]*S2[a] - gamma3*I3[a] - (um+u[a])*I3[a];
       dydt[a+7*agegroups] =  gamma3*I3[a] + gamma3*I4[a] - sigma3*lambda[a]*S3[a] - (um+u[a])*S3[a];
       dydt[a+8*agegroups] =  sigma3*lambda[a]*S3[a] - gamma3*I4[a] - (um+u[a])*I4[a];
       
       if (a >1 ){
		    dydt[a] = dydt[a]+ u[a-1]*M[a-1];
		    dydt[a+agegroups] = dydt[a+agegroups] + u[a-1]*S0[a-1];
		    dydt[a+2*agegroups] = dydt[a+2*agegroups] + u[a-1]*I1[a-1];
		    dydt[a+3*agegroups] =dydt[a+3*agegroups]+u[a-1]*S1[a-1];
		    dydt[a+4*agegroups] = dydt[a+4*agegroups]+u[a-1]*I2[a-1];
		    dydt[a+5*agegroups] = dydt[a+5*agegroups]+u[a-1]*S2[a-1];
		    dydt[a+6*agegroups] = dydt[a+6*agegroups]+u[a-1]*I3[a-1];
		    dydt[a+7*agegroups] = dydt[a+7*agegroups]+u[a-1]*S3[a-1];
		    dydt[a+8*agegroups] = dydt[a+8*agegroups]+u[a-1]*I4[a-1];
		}
       
       }
       
       
      
      return dydt;
  }
}

// Fixed data is declared in the data block:
data {
  int<lower=1> n_days;
  int<lower=1> agegroups;
  int q;
  real y0[9*agegroups];
  real t0;
  real ts[n_days];
  int N;
  real birthrate[agegroups];
  real um;
  real rho1;
  real rho2;
  real gamma1;
  real gamma2;
  real gamma3;
  real sigma1;
  real sigma2;
  real sigma3;
  int<lower=1> hosp_cases[n_days];
  real hosp1;
  real hosp2;
  real hosp3;
  real omega;
  real c2[agegroups*agegroups];
  real u[agegroups];
}


transformed data {
  real x_r[12];
  int x_i[2];
  x_i[1] = agegroups;
  x_i[2] = q;
  
  x_r[1:agegroups]= birthrate;
  x_r[agegroups+1]=um;
  x_r[agegroups+2]=rho1;
  x_r[agegroups+3]=rho2;
  x_r[agegroups+4]=gamma1;
  x_r[agegroups+5]=gamma2;
  x_r[agegroups+6]=gamma3;
  x_r[agegroups+7]=sigma1;
  x_r[agegroups+8]=sigma2;
  x_r[agegroups+9]=sigma3;
  x_r[agegroups+10]=omega;
  x_r[(agegroups+11):(11+2*agegroups)]=u;
  x_r[(12+2*agegroups):(12+2*agegroups+agegroups*agegroups)]=c2;
}

parameters {
  real<lower=0> beta;
  real<lower=0,upper=1> b1;
  real trans;
}

transformed parameters{
  real y[n_days, 9*agegroups]; // y is the states matrix that has row length=n_days and columns=9.
  real<lower=0,upper=2*pi()> phi=(2*pi()*exp(trans))/(1+exp(trans));
  real rel_inf[n_days];
  real season_beta[n_days];
    // outcomes
  vector<lower = 0>[n_days] output_hosp; // overall case incidence by day
  {
    real theta[3];
    theta[1] = beta;
    theta[2] = b1;
    theta[3] = trans;
    
    y = integrate_ode_rk45(msis, y0, t0, ts, theta, x_r, x_i);

  }
  
  for(i in 1:n_days){
    
    for(j in 1:9){
    if (y[i,j] <= 0.0) y[i,j] = 1e-12;}
    
    rel_inf[i]=(y[i,3]+rho1*y[i,5]+rho2*y[i,7]+rho2*y[i,9])/sum(y[i,]);
    
    season_beta[i]=(1+b1*cos(2*pi()*i/12 - phi))*beta;
    
    output_hosp[i] = hosp1*y[i,2]*rel_inf[i]*season_beta[i]+hosp2*sigma1*y[i,4]*rel_inf[i]*season_beta[i] +hosp3*sigma2*y[i,6]*rel_inf[i]*season_beta[i] +hosp3*sigma3*y[i,8]*rel_inf[i]*season_beta[i];
  }
}

model {
  //priors
  beta ~ normal(13.7,10); //truncated at 0
  b1 ~ normal(0.52, 0.25); //truncated at 0
  trans ~ normal(0.2, 1); //truncated at 0
  
  for(i in 1:n_days) {
  //sampling distribution
  target += poisson_lpmf(hosp_cases[i] | output_hosp[i]);
  }
}

