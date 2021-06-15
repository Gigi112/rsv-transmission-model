function [Base_site, Routine_site, RoutineC05_site, RoutineC15_site, RoutineC25_site, RoutineC99_site, ParamsDis_site, ParamsVax_site] = fn_vax_sim(output, vax_samp, agepar, params, randout, randout2, population, int_plans, simulationnames, sites, draw, i)    
    
    popdist = population.(sites{i})./sum(population.(sites{i}));
    pop100k = population.(sites{i})./sum(population.(sites{i}))*1e5;
    
for j=1:draw
        % Call the vaccination parameter. Makes sure all sites get the same vax_pars.
        randnow2 = randout2(j);
        vacpar.omega_v = -log(1-vax_samp.omega(randnow2))/52; % turn yearly to weekly rate
        vacpar.veff = vax_samp.nu(randnow2);

        randnow = randout(j);
        estpar.beta = output.(sites{i}).beta(randnow);
        estpar.mult1 = exp(-output.(sites{i}).logm1(randnow));
        estpar.mult2 = 1-exp(-output.(sites{i}).logm2(randnow));
        estpar.prop = output.(sites{i}).prop(randnow);
        estpar.r = 1; % r.(sites{i}{i})(randi);
        estpar.rC = output.(sites{i}).r(randnow);
        estpar.sensrate = output.(sites{i}).sensrate(randnow);
        params.epsilon = 0;
        % R0 = beta/(mu+delta)*(1+r*theta*delta/mu);
        params.R0 = estpar.beta/(agepar.(sites{i}).mu'*popdist+params.delta)*(1+(estpar.rC*(agepar.(sites{i}).theta'*popdist)*params.delta)./(agepar.(sites{i}).mu'*popdist));

        ParamsDis_site(j) = estpar;        
        ParamsVax_site(j).omegav = vacpar.omega_v;        
        ParamsVax_site(j).veff = vacpar.veff;        

        for k=1:6
        % include parameters about the age of vaccination...
        vacpar.v1 = zeros(1+30*52,11);
        vacpar.v1(5:end, int_plans(k,1)) = int_plans(k,3);
        vacpar.massvacc = zeros(1+30*52,11);
        vacpar.massvacc(1:4, int_plans(k,1):int_plans(k,2)) = min((1-(1-int_plans(k,4))^.25),.62); 
      
       % Simulation equation
       if i==4
        tspan=200*52+1;
       else
        tspan=100*52+1;
       end
       
        [out, pop17] = fn_forward_sim_vacc(params, agepar.(sites{i}), estpar, vacpar, pop100k, tspan, sites{i}); 

        eval([simulationnames{k} '(j) = out;']);        
        
        end
    end
end
    