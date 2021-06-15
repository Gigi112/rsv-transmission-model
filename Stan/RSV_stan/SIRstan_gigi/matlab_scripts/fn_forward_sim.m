function [eqmcases, eqmchronic, eqmpop] = fn_forward_sim(params, agepar, estpar, population, tspan, site)
% Forward simulations of typhoid system in Kolkata
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

mub = agepar.mub;
u = agepar.u;
mu = agepar.mu;
theta = agepar.theta;
theta2 = agepar.theta2;

omega = params.omega;
% r = params.r;
% rC = params.rC;
alpha = params.alpha;
delta = params.delta;
epsilon = params.epsilon;

% Population in year -50
pop0 = [population; zeros(al*6,1)];
% throw an infected in one of the four categories to start the simulation
ageptzero = 3; %randsample(1:al, 1, true, pop0(1:al)/sum(pop0(1:al)));
pop0(ageptzero) = pop0(ageptzero)-1;
pop0(al+ageptzero) = pop0(al+ageptzero)+1;

% Simulate
% The reason that the simulation is not more straight forward is that we
% want to keep the population size and age-structure similar to the one we
% obtained from the literature. If we were to run the model to equilibrium, 
% the age structure of the population would shift significantly (as is 
% happening in the real world). Therefore we run the model for 5 years, and 
% then assign the proportion of the population in each of the disease 
% compartments into a population with an age-structure equal to the 
% age-structure that we observe.

options = odeset('RelTol',1e-6);
popi = pop0;

cycles = tspan/tspan;
fitp = zeros(cycles, al*6);
fitc = zeros(cycles, al);
fitch = zeros(cycles, al);

try
for tee = 1:cycles
    [t, pop] = ode45(@(t, pop) fn_SIR(t, pop, betap, mub, u, mu, omega, r, rC, alpha, delta, theta, theta2, epsilon, al, population), 1:tspan, popi, options); % kappa,
    fitp(tee, 1:(al*6)) = pop(end, 1:(al*6));

    for this = 1:al
        seg = this:al:(5*al+this);
        popi(seg)= pop0(this)*fitp(tee,seg)/sum(fitp(tee,seg));
            % How much of 'this' age compartment is is each particular compartment (S, Im, R, Sq, Iq, C)
    end
    fitc(tee,:) = pop(end, (6*al+1):(7*al)) - pop(end-52, (6*al+1):(7*al));
    fitch(tee,:) = pop(end, (5*al+1):(6*al));
end
catch
    sprintf('ODE error: ', site)
    fitc = zeros(cycles, al);
    fitch = zeros(cycles, al);    
    fitp = zeros(cycles, al*6);
end

eqmchronic = fitch(cycles,:)';
eqmcases = fitc(cycles,:)';
eqmpop = ones(al, 1); % what's my pop distribution at the end? 

for this = 1:al
    seg = this:al:(5*al+this);
    eqmpop(this, 1) = sum(fitp(cycles, seg));
end
    
end
