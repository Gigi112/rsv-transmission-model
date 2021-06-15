function [out] = fn_vax_sim_kolkata(ParamsDis, ParamsVax, prices, discountlevel, agepar, params, population, Sim11ages, z, draw)    
    
for j=1:draw
                
        betanew = ParamsDis(j,1).beta;
        mult1new = ParamsDis(j,1).mult1;
        mult2new = ParamsDis(j,1).mult2;
        r = 1;
        rC = ParamsDis(j,1).rC;

        al = length(population);

        betap = betanew*ones(al, 1);
        betap(agepar.under5==1) = mult1new*betanew;
        betap(agepar.under5==2) = mult2new*betanew;

        mub = agepar.mub;
        u = agepar.u;
        mu = agepar.mu;
        mu_real = agepar.mu_real;
        theta = agepar.theta;
        theta2 = agepar.theta2;

        omega = params.omega;
        alpha = params.alpha;
        delta = params.delta;
        epsilon = params.epsilon;

        veff = ParamsVax(j,1).veff; 
        omega_v = ParamsVax(j,1).omega_v; 
        
        % Population in year 0 from Sim11ages
        pop17 = [Sim11ages{j,1}(end, 1:al*6)./sum(Sim11ages{j,1}(end, 1:al*6))*1e5, zeros(1,al*12)];
        
        % time vector for ode
        tspan = 31*52+4+1;
        % time vector for cost-effectiveness:
        obs_time = 5:52:(52*31+4+1);
        
        options = odeset('RelTol',1e-6);

    for pricelevel=1:length(prices)

        % 3 price levels: $1,$2,$5,$105.88. Add administration, 
        % supply fees: $3.55 + $0.17 = 3.72. AND 15*pricelevel wastage. 
        % 5 discounts: 0%, 25%, 50%, 75%, 100% 
        % 1 vs 2 doses (for children <5; <5 campaign is free anyways?)

      for doses=1:2
      totalprice = (1.15*prices(pricelevel)+3.55+0.17)*doses;
      
      for discount=1:length(discountlevel)
      finalprice=totalprice*(1-discountlevel(discount));
          
      if discountlevel(discount)==1
      vcov1 = ParamsVax(j,1).coverageatfree;
      else
      vcov1 = ParamsVax(j,1).coverageatfree*(exp(-log10(finalprice)/ParamsVax(j,1).scale)^ParamsVax(j,1).shape); 
      end
      
      vcov2 = ParamsVax(j,1).coverageatfree;

    int_plans = [1, 1, 0, 0; ... % No intervention
        2, 2, vcov1, 0; ... % Routine vaccination (the 3 in the second position doesn't matter)
        2, 3, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 5 years of age
        2, 5, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 15 years of age
        4, 4, vcov1, 0; ... % Routine vaccination at 5 years of age
        4, 5, vcov1, vcov2]; % Routine vaccination at 5 years of age and catch up campaign for people up to 15 years of age

    v1 = zeros(tspan,11); % when routine vaccination should take place (time & age)
    v1(5:end, int_plans(z,1)) = int_plans(z,3);

    massvacc = zeros(tspan,11); % when mass vaccination should take place (time & age)
    massvacc(1:4,  int_plans(z,1):int_plans(z,2)) = min(1-(1-int_plans(z,4))^0.25, 0.62);
    
    [t,y] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata(t, pop, betap, mub, u, mu_real, omega, omega_v, r, rC, alpha, delta, theta, theta2, epsilon, v1, veff, massvacc, al, population), 1:tspan, pop17, options);

    % POPULATION in each age group each year we are observing them...
    out{j,pricelevel}{doses, discount}.pop = ones(length(obs_time), al);
    out{j,pricelevel}{doses, discount}.popu = ones(length(obs_time), al);
    out{j,pricelevel}{doses, discount}.popv = ones(length(obs_time), al);

    for this = 1:al
     seg = [this:al:(13*al+this)] ;
     out{j,pricelevel}{doses, discount}.pop(:, this) = sum(y(obs_time, seg),2);

     seg = [this:al:(5*al+this)] ;
     out{j,pricelevel}{doses, discount}.popu(:, this) = sum(y(obs_time, seg),2);

     seg = [(6*al):al:(13*al+this)] ;
     out{j,pricelevel}{doses, discount}.popv(:, this) = sum(y(obs_time, seg),2);
    end

    % INCIDENCE AMONG PEOPLE WHO WERE NOT VACCINATED.
             out{j,pricelevel}{doses, discount}.cumI1 = y(obs_time, (16*al+1):(17*al));
    % INCIDENCE AMONG VACCINATED PEOPLE
             out{j,pricelevel}{doses, discount}.cumI1v = y(obs_time, (17*al+1):(18*al));
    % PREVALENT CASES OF CHRONIC CARRIERS
             out{j,pricelevel}{doses, discount}.chr = y(obs_time, (5*al+1):(6*al));
             out{j,pricelevel}{doses, discount}.chrv = y(obs_time, (13*al+1):(14*al));
    % VACCINES ADMINISTERED (CUMULATIVE VACCINES)
             out{j,pricelevel}{doses, discount}.cumdosesr = y(obs_time, (14*al+1):(15*al));
    % VACCINES ADMINISTERED FOR CATCH-UP & MASS CAMPAIGNS (CUMULATIVE VACCINES)
             out{j,pricelevel}{doses, discount}.cumdosesc = y(obs_time, (15*al+1):(16*al));
    % PEOPLE EFFECTIVELY VACCINATED at any point in time
             out{j,pricelevel}{doses, discount}.V1 = y(obs_time, 6*al+1:7*al) ; %V1 + V2
             out{j,pricelevel}{doses, discount}.V2 = y(obs_time, 7*al+1:8*al) ; %V1 + V2
             out{j,pricelevel}{doses, discount}.int_plans = int_plans(z,:);
      end
      end
    end
    
    end
end
    