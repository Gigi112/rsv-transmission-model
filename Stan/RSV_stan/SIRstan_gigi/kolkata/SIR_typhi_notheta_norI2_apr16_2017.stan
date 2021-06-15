functions {
       real[] make_x_r(int agegroups,
                real population, 
                real mub,
                // real r,
                real alpha,
                real omega,
                real delta,
                real[] mu,
                real[] theta,
                real[] u) {
                real x_r[agegroups*8];
                    x_r[1] = population; // total population size
                    x_r[1+agegroups] = mub; // additions due to births only to the youngest age group
                    for (j in 2:agegroups) {
                    x_r[j] = 0; // total population size is not age dependent...
                    x_r[j+agegroups] = 0; // additions due to births are zero in all but the youngest group
                    }
                    for (i in 1:agegroups) {
                    x_r[i+agegroups*2] = alpha; // case fatality rate of typhi
                    x_r[i+agegroups*3] = omega; // rate of waning immunity
                    x_r[i+agegroups*4] = delta; // rate at which an infectious individual stops transmitting typhi
                    // x_r[i+agegroups*8] = r_equality; // relative contribution of asymptomatic carriers
                    }
                    x_r[(1+agegroups*5):(agegroups*6)] = mu; // death rate, by age
                    x_r[(1+agegroups*6):(agegroups*7)] = theta; // probability of becoming a chronic carrier, by age
                    x_r[(1+agegroups*7):(agegroups*8)] = u; // rate of aging out of each age group, by age
                return x_r;
                    } 
       real[] make_params(real R0, 
                real logm1,
                real logm2,
                real r,
                int agegroups,
		 real mub,
		 real delta,
		 real population) {
                real params[agegroups+1];
                params[1] = exp(-logm1)*R0*(mub/population+delta)/(1+delta*r*0.01/mub/population); 
                params[2] = (1-exp(-logm2))*R0*(mub/population+delta)/(1+delta*r*0.01/mub/population);
                for (i in 3:agegroups)
                params[i] = R0*(mub/population+delta)/(1+delta*r*0.01/mub/population); // 0.01 is for average theta throughout population.
                params[agegroups+1] = r;
               return params;
                }
       real[] sir(real t,
           real[] y, // these will be the initial population distribution
           real[] theta, // these will be the parameters to be estimated (except reporting and blood sensitivity)
           real[] x_r, // x_r will be the fixed parameters
           int[] x_i) { // x_i will be the number of age groups, whether I2 has zero chance of becoming a chronic carrier, and whether I2 contributes the same as I1 or as C.
            real dydt[8*x_i[1]]; 
            for (j in 1:1) {
                dydt[j] = x_r[j+x_i[1]] - y[j+x_i[1]*7]*y[j] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j]; 
					// primary susceptibles: births - infections - (aging out + background deaths)
                dydt[j+x_i[1]] = y[j+x_i[1]*7]*y[j] - x_r[j+x_i[1]*4]*y[j+x_i[1]] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]]; 
					// first-time infections: infections - recovery - (aging out + background deaths)
                dydt[j+x_i[1]*2] = x_r[j+x_i[1]*4]*(1-x_r[j+x_i[1]*2]-x_r[j+x_i[1]*6])*y[j+x_i[1]] + x_r[j+x_i[1]*4]*(1-x_i[2]*x_r[j+x_i[1]*6])*y[j+x_i[1]*4] - x_r[j+x_i[1]*3]*y[j+x_i[1]*2] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*2]; // recovered
                dydt[j+x_i[1]*3] = x_r[j+x_i[1]*3]*y[j+x_i[1]*2] - y[j+x_i[1]*7]*y[j+x_i[1]*3] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*3]; // secondary susceptibles
                dydt[j+x_i[1]*4] = y[j+x_i[1]*7]*y[j+x_i[1]*3] - x_r[j+x_i[1]*4]*y[j+x_i[1]*4] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*4]; // secondary infectious
                dydt[j+x_i[1]*5] = x_r[j+x_i[1]*4]*x_r[j+x_i[1]*6]*(y[j+x_i[1]] + x_i[2]*y[j+x_i[1]*4]) - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*5]; // chronic cases
                dydt[j+x_i[1]*6] = y[j+x_i[1]*7]*y[j]; // cumulative incidence: cumulative sum of first-time infections
                dydt[j+x_i[1]*7] = theta[j]*(sum(y[(x_i[1]+1):(2*x_i[1])])+(theta[x_i[1]+1]*x_i[3]+(1-x_i[3]))*sum(y[(4*x_i[1]+1):(5*x_i[1])])+theta[x_i[1]+1]*sum(y[(5*x_i[1]+1):(6*x_i[1])]))/x_r[1] - y[j+x_i[1]*7]; 
                        // LAMBDA: beta*[(first-time infections) + r*(indicator for 1 or fraction)*(second-time infections) + r*chronics]/population
            }
            for (j in 2:x_i[1]) { // these compartments include aging IN
                dydt[j] = - y[j+x_i[1]*7]*y[j] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j] + x_r[j-1+x_i[1]*7]*y[j-1]; // primary susceptibles
                dydt[j+x_i[1]] = y[j+x_i[1]*7]*y[j] - x_r[j+x_i[1]*4]*y[j+x_i[1]] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]] + x_r[j-1+x_i[1]*7]*y[j-1+x_i[1]]; // primary infections
                dydt[j+x_i[1]*2] = x_r[j+x_i[1]*4]*(1-x_r[j+x_i[1]*2]-x_r[j+x_i[1]*6])*y[j+x_i[1]] + x_r[j+x_i[1]*4]*(1-x_i[2]*x_r[j+x_i[1]*6])*y[j+x_i[1]*4] - x_r[j+x_i[1]*3]*y[j+x_i[1]*2] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*2] + x_r[j-1+x_i[1]*7]*y[j-1+x_i[1]*2]; // recovered
                dydt[j+x_i[1]*3] = x_r[j+x_i[1]*3]*y[j+x_i[1]*2] - y[j+x_i[1]*7]*y[j+x_i[1]*3] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*3] + x_r[j-1+x_i[1]*7]*y[j-1+x_i[1]*3]; // secondary susceptibles
                dydt[j+x_i[1]*4] = y[j+x_i[1]*7]*y[j+x_i[1]*3] - x_r[j+x_i[1]*4]*y[j+x_i[1]*4] - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*4] + x_r[j-1+x_i[1]*7]*y[j-1+x_i[1]*4]; // secondary infectious
                dydt[j+x_i[1]*5] = x_r[j+x_i[1]*4]*x_r[j+x_i[1]*6]*(y[j+x_i[1]] + x_i[2]*y[j+x_i[1]*4]) - (x_r[j+x_i[1]*7]+x_r[j+x_i[1]*5])*y[j+x_i[1]*5] + x_r[j-1+x_i[1]*7]*y[j-1+x_i[1]*5]; // chronic cases
                dydt[j+x_i[1]*6] = y[j+x_i[1]*7]*y[j]; // cumulative incidence
                dydt[j+x_i[1]*7] = theta[j]*(sum(y[(x_i[1]+1):(2*x_i[1])])+(theta[x_i[1]+1]*x_i[3]+(1-x_i[3]))*sum(y[(4*x_i[1]+1):(5*x_i[1])])+theta[x_i[1]+1]*sum(y[(5*x_i[1]+1):(6*x_i[1])]))/x_r[1] - y[j+x_i[1]*7]; // LAMBDA
                }
           return dydt;
           }
    }

data { 
    int<lower=1> ts; 
    real<lower=0> tspan[ts];
    int<lower=0> agegroups; // Number of age groups in each site
	int<lower=0, upper=1> thetaindicator;
	int<lower=0, upper=1> rindicator;
    int<lower=0> cases[agegroups];
    int<lower=0> casesdhaka2[2];
    real alpha;
    real omega; 
    real delta; 
    real population; 
    real mub;
    int<lower=0, upper=1> dhaka;

    // Age specific fixed inputs and parameters.
    real<lower=0> y0[agegroups*8];
    int<lower=0, upper=agegroups*6> y0index[agegroups, 6];
    real<lower=-1, upper=1> mu[agegroups];
    real<lower=0> u[agegroups];
    real<lower=0, upper=1> theta[agegroups];
    real<lower=0, upper=1> participation[agegroups];
    real<lower=0> vol[agegroups];
}

transformed data {
   real t0; 
   real x_r[agegroups*8];
   int x_i[3];

   x_i[1] = agegroups;
   x_i[2] = thetaindicator;
   x_i[3] = rindicator;

   t0 = 0;
   x_r = make_x_r(agegroups, population, mub, alpha, omega, delta, mu, theta, u);
}

parameters { 
    real<lower=1, upper=2.5> R0;
    real<lower=0> logm1;
    real<lower=0> logm2;
    real<lower=0> sensrate;
    real<lower=0> logprop;
    real<lower=0, upper=1> r; 
}

transformed parameters {  
       real params[agegroups+1];
       real blood[agegroups];
       real prop;

       params = make_params(R0, logm1, logm2, r, agegroups, mub, delta, population);


       for (i in 1:agegroups) 
            blood[i] = 1-exp(-sensrate*vol[i]);

       prop = 1-exp(-logprop);

       }

model {
       real yhat[ts, agegroups*8];
       real yhatpop[agegroups];
    // Disease process parameters
        // Local parameters
        R0 ~ uniform(1, 3.5);

        logm1 ~ exponential(1); // uniform(0, 1);
        logm2 ~ exponential(1); // uniform(0, 1);
        r ~ gamma(6.7194, 1/0.0517); // beta(6.34, 19.4); // Same parameterization as in matlab. 
    
    // Run ODE
       yhat = integrate_ode_rk45(sir, y0[1:agegroups*8], t0, tspan, params, x_r, x_i);
    for (i in 1:agegroups)
        yhatpop[i] = sum(yhat[ts, y0index[i,1:6]]);

    // Observation process parameters
        // Local parameters
        logprop ~ exponential(1);
        sensrate ~ gamma(44.8086, 1/0.0031); // gamma's beta in stan is parameterized as the inverse of the gamma's beta in matlab

    // Likelihood
    for (i in 1:agegroups)
        cases[i] ~ poisson((yhat[ts, (i+agegroups*6)] - yhat[ts-52, (i+agegroups*6)])/yhatpop[i]*y0[i]*participation[i]*blood[i]*prop); 
    if (dhaka==1){
        casesdhaka2[1] ~ poisson((sum(yhat[ts, (1+agegroups*6):(2+agegroups*6)]) - sum(yhat[ts-52, (1+agegroups*6):(2+agegroups*6)]))*participation[2]*blood[2]*prop);
        casesdhaka2[2] ~ poisson((sum(yhat[ts, (3+agegroups*6):(4+agegroups*6)]) - sum(yhat[ts-52, (3+agegroups*6):(4+agegroups*6)]))*participation[4]*blood[4]*prop);
    }

} // generated quantities (for WAIC, R0, etc) might be a good idea later.
