function [Base_site, ParamsDis_site] = fn_vax_sim_nointervention(output, agepar, params, randout, population, draw, tspan)    

    popdist = population./sum(population);
    pop100k = population./sum(population)*1e5;
    
for j=1:draw
        randnow = randout(j);
        estpar.beta = output.beta(randnow);
        estpar.mult1 = exp(-output.logm1(randnow));
        estpar.mult2 = 1-exp(-output.logm2(randnow));
        estpar.prop = output.prop(randnow);
        estpar.r = 1; % r.(sites{i}{i})(randi);
        estpar.rC = output.r(randnow);
        estpar.sensrate = output.sensrate(randnow);
        params.epsilon = 0;
        % R0 = beta/(mu+delta)*(1+r*theta*delta/mu);
        params.R0 = estpar.beta/(agepar.mub(1)+params.delta)*(1+(estpar.rC*(agepar.theta'*popdist)*params.delta))./agepar.mub(1);
        
        ParamsDis_site(j) = estpar;        
        
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
        alpha = params.alpha ;
        delta = params.delta ;
        epsilon = params.epsilon ;

        % Population in year 0, eqmpop has zeros in vaccination compartments

        pop0 = [population; zeros(al*6,1)];
        % throw an infected in one of the four categories to start the simulation
        ageptzero = 3; %randsample(1:al, 1, true, pop0(1:al)/sum(pop0(1:al)));
        pop0(ageptzero) = pop0(ageptzero)-1;
        pop0(al+ageptzero) = pop0(al+ageptzero)+1;

        options = odeset('RelTol',1e-6);
         
        [t, pop] = ode45(@(t, pop) fn_SIR(t, pop, betap, mub, u, mu, omega, r, rC, alpha, delta, theta, theta2, epsilon, al, population), 1:tspan, pop0, options); % kappa,
         
       % Simulation equation      
        Base_site{j,1} = pop;        
        
    end
end
    