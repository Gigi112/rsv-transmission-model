function [out, pop17] = fn_forward_sim_vacc(params, agepar, estpar, vacpar, population, tspan, site)

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
v1 = vacpar.v1; % when routine vaccination should take place (time & age)
massvacc = vacpar.massvacc; % when mass vaccination should take place (time & age)
omega_v = vacpar.omega_v; 
    
% Population in year 0, eqmpop has zeros in vaccination compartments

pop0 = [population; zeros(al*6,1)];
% throw an infected in one of the four categories to start the simulation
ageptzero = 3; %randsample(1:al, 1, true, pop0(1:al)/sum(pop0(1:al)));
pop0(ageptzero) = pop0(ageptzero)-1;
pop0(al+ageptzero) = pop0(al+ageptzero)+1;

options = odeset('RelTol',1e-6);
         
[t, pop] = ode45(@(t, pop) fn_SIR(t, pop, betap, mub, u, mu, omega, r, rC, alpha, delta, theta, theta2, epsilon, al, population), 1:tspan, pop0, options); % kappa,
         fitp(1, 1:(al*6)) = pop(end, 1:(al*6));
         
for this = 1:al
    seg = this:al:(5*al+this);
    popi(seg)= pop0(this)*fitp(1,seg)/sum(fitp(1,seg));
    % How much of 'this' age compartment is is each particular compartment (S, Im, R, Sq, Iq, C)
end
pop17 = [popi(1:(6*al))'; zeros(al*12,1)];

[t,y] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata(t, pop, betap, mub, u, mu_real, omega, omega_v, r, rC, alpha, delta, theta, theta2, epsilon, v1, veff, massvacc, al, population), 1:(1+52*30), pop17, options);
 
% time vector for cost-effectiveness:
obs_time = 1:52:(52*30+1);
         
% POPULATION in each age group each year we are observing them...
out.pop = ones(length(obs_time), al);
out.popu = ones(length(obs_time), al);
out.popv = ones(length(obs_time), al);

for this = 1:al
         seg = [this:al:(13*al+this)] ;
         out.pop(:, this) = sum(y(obs_time, seg),2);
         
         seg = [this:al:(5*al+this)] ;
         out.popu(:, this) = sum(y(obs_time, seg),2);
         
         seg = [(6*al):al:(13*al+this)] ;
         out.popv(:, this) = sum(y(obs_time, seg),2);
end

% INCIDENCE AMONG PEOPLE WHO WERE NOT VACCINATED.
         out.cumI1 = y(obs_time, (16*al+1):(17*al));
% INCIDENCE AMONG VACCINATED PEOPLE
         out.cumI1v = y(obs_time, (17*al+1):(18*al));
% PREVALENT CASES OF CHRONIC CARRIERS
         out.chr = y(obs_time, (5*al+1):(6*al));
         out.chrv = y(obs_time, (13*al+1):(14*al));
% VACCINES ADMINISTERED (CUMULATIVE VACCINES)
         out.cumdosesr = y(obs_time, (14*al+1):(15*al));
% VACCINES ADMINISTERED FOR CATCH-UP & MASS CAMPAIGNS (CUMULATIVE VACCINES)
         out.cumdosesc = y(obs_time, (15*al+1):(16*al));
% PEOPLE EFFECTIVELY VACCINATED at any point in time
         out.V1 = y(obs_time, 6*al+1:7*al) ; %V1 + V2
         out.V2 = y(obs_time, 7*al+1:8*al) ; %V1 + V2

end
