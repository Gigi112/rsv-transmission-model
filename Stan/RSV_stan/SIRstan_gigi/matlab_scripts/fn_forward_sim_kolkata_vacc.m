function out = fn_forward_sim_kolkata_vacc(params, agepar, estpar, population, tspan)

fn_vax_sim_kolkata(output, vax_samp, agepar, params, randout, randout2, population, int_plans, simulationnames, sites, draw, i);

% Produces a vector of the symptomatic cases at the end of one year at
% equillibrum (incidence) and a vector of the sizes of each compartment 
% at the end of a year in equillibrium (prevalence in each compartment). 

betanew = estpar.beta;
mult1new = estpar.mult1;
mult2new = estpar.mult2;
r = estpar.r;
rC = estpar.rC;

al = length(population);

betap = betanew*ones(al, 1);
betap(agepar.under5==1) = mult1new*betanew ;
betap(agepar.under5==2) = mult2new*betanew ;

mub = agepar.mub ;
u = agepar.u ;
mu = agepar.mu ;
mu_real = agepar.mu_real ;
theta = agepar.theta ;
theta2 = agepar.theta2 ;

omega = params.omega ;
% r = params.r ;
% rC = params.rC ;
alpha = params.alpha ;
delta = params.delta ;
epsilon = params.epsilon ;

veff = vacpar.veff; 
omega_v = vacpar.omega_v; 
  
% Population in year 0 from Sim11ages

% time vector for cost-effectiveness:
obs_time = 1:52:(52*30+1);

% Cycle through all 6 prices/coverage combinations
% Unless this is BASE

for j=1:length(vacpar.prices)

[t,y] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata(t, pop, betap, mub, u, mu_real, omega, omega_v, r, rC, alpha, delta, theta, theta2, epsilon, v1, veff, massvacc, al, population), 1:(1+52*30), pop17, options);

v1 = vacpar.v1; % when routine vaccination should take place (time & age)
massvacc = vacpar.massvacc; % when mass vaccination should take place (time & age)
  
% POPULATION in each age group each year we are observing them...
out{j,1}.pop = ones(length(obs_time), al);
out{j,1}.popu = ones(length(obs_time), al);
out{j,1}.popv = ones(length(obs_time), al);

for this = 1:al
         seg = [this:al:(13*al+this)] ;
         out{j,1}.pop(:, this) = sum(y(obs_time, seg),2);
         
         seg = [this:al:(5*al+this)] ;
         out{j,1}.popu(:, this) = sum(y(obs_time, seg),2);
         
         seg = [(6*al):al:(13*al+this)] ;
         out{j,1}.popv(:, this) = sum(y(obs_time, seg),2);
end

% INCIDENCE AMONG PEOPLE WHO WERE NOT VACCINATED.
         out{j,1}.cumI1 = y(obs_time, (16*al+1):(17*al));
% INCIDENCE AMONG VACCINATED PEOPLE
         out{j,1}.cumI1v = y(obs_time, (17*al+1):(18*al));
% PREVALENT CASES OF CHRONIC CARRIERS
         out{j,1}.chr = y(obs_time, (5*al+1):(6*al));
         out{j,1}.chrv = y(obs_time, (13*al+1):(14*al));
% VACCINES ADMINISTERED (CUMULATIVE VACCINES)
         out{j,1}.cumdosesr = y(obs_time, (14*al+1):(15*al));
% VACCINES ADMINISTERED FOR CATCH-UP & MASS CAMPAIGNS (CUMULATIVE VACCINES)
         out{j,1}.cumdosesc = y(obs_time, (15*al+1):(16*al));
% PEOPLE EFFECTIVELY VACCINATED at any point in time
         out{j,1}.V1 = y(obs_time, 6*al+1:7*al) ; %V1 + V2
         out{j,1}.V2 = y(obs_time, 7*al+1:8*al) ; %V1 + V2
end

end
