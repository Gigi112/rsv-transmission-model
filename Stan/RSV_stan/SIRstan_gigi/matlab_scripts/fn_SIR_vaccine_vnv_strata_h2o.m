function dydt = fn_SIR_vaccine_vnv_strata_h2o(t, pop, betap, mub, u, mu, omega, omega_v, r, rC, alpha, delta, theta, theta2, epsilon, v1, veff, massvacc, al, population)
%Differential equations for typhoid model

% make two more parameters, and then we are ready for ode's
lambdap=betap(:,round(t))*sum(pop((al+1):(2*al))+pop((9*al+1):(10*al))+r*(pop((4*al+1):(5*al))+pop((12*al+1):(13*al)))+ rC*(pop((5*al+1):(6*al))+ pop((13*al+1):(14*al))))/sum(population);
% somehow exceeds matrix dimensions
births = mub*sum(pop(1:(14*al)));

for i=1:al
% DISEASE PROCESS
    % disease processes for individuals who have NOT been vaccinated, ever.
    dydt(i,1) = births(i) + epsilon*pop(i+3*al) - lambdap(i)*pop(i) - (u(i)+mu(i))*pop(i); %dS1/dt
    dydt(i+al,1) = lambdap(i)*pop(i) - delta*pop(i+al) - (u(i)+mu(i))*pop(i+al); %dI1/dt 
    dydt(i+2*al,1) = delta*((1-alpha-theta(i))*pop(i+al) + (1-theta2(i))*pop(i+4*al)) - omega*pop(i+2*al) - (u(i)+mu(i))*pop(i+2*al); %dR/dt
    dydt(i+3*al,1) = omega*pop(i+2*al) - epsilon*pop(i+3*al) - lambdap(i)*pop(i+3*al) - (u(i)+mu(i))*pop(i+3*al); %dS2/dt 
    dydt(i+4*al,1) = lambdap(i)*pop(i+3*al) - delta*pop(i+4*al) - (u(i)+mu(i))*pop(i+4*al); %dI2/dt
    dydt(i+5*al,1) = delta*(theta(i)*pop(i+al) + theta2(i)*pop(i+4*al)) - (u(i)+mu(i))*pop(i+5*al); %dC/dt
 
    % vaccinated states
    dydt(i+6*al,1) = -omega_v*pop(i+6*al) - (u(i)+mu(i))*pop(i+6*al); %dV1/dt
    dydt(i+7*al,1) = -omega_v*pop(i+7*al) - (u(i)+mu(i))*pop(i+7*al); %dV2/dt
    
    % disease process for individuals who have been vaccinated before.
    dydt(i+8*al,1) = omega_v*pop(i+6*al) + epsilon*pop(i+11*al) - lambdap(i)*pop(i+8*al) - (u(i)+mu(i))*pop(i+8*al); %dS1v/dt
    dydt(i+9*al,1) = lambdap(i)*pop(i+8*al) - delta*pop(i+9*al) - (u(i)+mu(i))*pop(i+9*al); %dI1v/dt 
    dydt(i+10*al,1) = delta*((1-alpha-theta(i))*pop(i+9*al) + (1-theta2(i))*pop(i+12*al)) - omega*pop(i+10*al) - (u(i)+mu(i))*pop(i+10*al); %dRv/dt
    dydt(i+11*al,1) = omega*pop(i+10*al) + omega_v*pop(i+7*al) - epsilon*pop(i+11*al) - lambdap(i)*pop(i+11*al) - (u(i)+mu(i))*pop(i+11*al); %dS2v/dt 
    dydt(i+12*al,1) = lambdap(i)*pop(i+11*al) - delta*pop(i+12*al) - (u(i)+mu(i))*pop(i+12*al); %dI2v/dt
    dydt(i+13*al,1) = delta*(theta(i)*pop(i+9*al) + theta2(i)*pop(i+12*al)) - (u(i)+mu(i))*pop(i+13*al); %dCv/dt
     
    % dydt(i+14*al,1) keeps track of vaccines distributed in routine vaccination campaigns.
    % dydt(i+15*al,1) keeps track of the vaccines distributed in mass and catch-up campaigns.
    
    % Cumulative incidence of symptomatic infection among unvaccinated individuals 
    dydt(i+16*al,1) = lambdap(i)*pop(i); % because I need the cumulative cases over a year, not just the cases in one point in time.
    % Cumulative incidence of infection among vaccinated individuals 
    dydt(i+17*al,1) = lambdap(i)*pop(i+8*al);

% VACCINES DURING ROUTINE IMMUNIZATION (AND THE AGING PROCESS)    
    if i>1 %aging + vaccination (no one vaccinated in age group=1
        % Compartments of individuals who have never had vaccination and 
        % are aging in. Some change to vaccinated comapartments at the time
        % of aging in.
        dydt(i,1) = dydt(i,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1); % dS1/dt
        dydt(i+al,1) = dydt(i+al,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1+al); % dI1/dt
        dydt(i+2*al,1) = dydt(i+2*al,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1+2*al); % dR/dt
        dydt(i+3*al,1) = dydt(i+3*al,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1+3*al); % dS2/dt
        dydt(i+4*al,1) = dydt(i+4*al,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1+4*al); % dI2/dt
        dydt(i+5*al,1) = dydt(i+5*al,1) + (1-v1(round(t),i))*u(i-1)*pop(i-1+5*al); % dC/dt
        
        % Vaccine compartments
        dydt(i+6*al,1) = dydt(i+6*al,1) + u(i-1)*pop(i-1+6*al) + veff*v1(round(t),i)*u(i-1)*(pop(i-1)); % dV1/dt    +pop(i-1+8*al) <- when we care about re-vaccinations
        dydt(i+7*al,1) = dydt(i+7*al,1) + u(i-1)*pop(i-1+7*al) + veff*v1(round(t),i)*u(i-1)*(pop(i-1+3*al)+pop(i-1+11*al)+pop(i-1+2*al)+pop(i-1+10*al)); % dV2/dt
        
        % Compartments designating aging of vaccinated individuals, and
        % entry of individuals newly vaccinated.
        dydt(i+8*al,1)=dydt(i+8*al,1) + u(i-1)*pop(i-1+8*al) + (1-veff)*v1(round(t),i)*u(i-1)*pop(i-1); % dS1v/dt
        dydt(i+9*al,1)=dydt(i+9*al,1) + u(i-1)*pop(i-1+9*al) + v1(round(t),i)*u(i-1)*pop(i-1+al); % dI1v/dt
        dydt(i+10*al,1)=dydt(i+10*al,1) + (1-veff*v1(round(t),i))*u(i-1)*pop(i-1+10*al) + (1-veff)*v1(round(t),i)*u(i-1)*pop(i-1+2*al); % dRv/dt
        dydt(i+11*al,1)=dydt(i+11*al,1) + (1-veff*v1(round(t),i))*u(i-1)*pop(i-1+11*al) + (1-veff)*v1(round(t),i)*u(i-1)*pop(i-1+3*al); % dS2v/dt
        dydt(i+12*al,1)=dydt(i+12*al,1) + u(i-1)*pop(i-1+12*al) + v1(round(t),i)*u(i-1)*pop(i-1+4*al); %dI2v/dt
        dydt(i+13*al,1)=dydt(i+13*al,1) + u(i-1)*pop(i-1+13*al) + v1(round(t),i)*u(i-1)*pop(i-1+5*al); %dCv/dt

        % Number of vaccines distributed during routine vaccination 
        dydt(i+14*al,1)=v1(round(t),i)*u(i-1)*sum(pop(i-1:al:(i-1+5*al))+pop((i-1+8*al):al:(i-1+13*al))); % dydt(i+14*al,1) + 
            % allows for revaccination for individuals who have come up to an age for booster doses.
            % Assumes no inefficiency? 
    end
    
% VACCINATIONS DURING ONE-TIME CAMPAIGNS: CATCH-UP AND MASS VACCINES
    % only valid for the 4 weeks before routine vaccination is to begin.
    if massvacc(round(t),i)>0 
        dydt(i,1)=dydt(i,1) - massvacc(round(t),i)*pop(i); %dS1/dt
        dydt(i+al,1)=dydt(i+al,1) - massvacc(round(t),i)*pop(i+al); %dI1/dt
        dydt(i+2*al,1)=dydt(i+2*al,1) - massvacc(round(t),i)*pop(i+2*al); %dR/dt
        dydt(i+3*al,1)=dydt(i+3*al,1) - massvacc(round(t),i)*pop(i+3*al); %dS2/dt
        dydt(i+4*al,1)=dydt(i+4*al,1) - massvacc(round(t),i)*pop(i+4*al); %dI2/dt
        dydt(i+5*al,1)=dydt(i+5*al,1) - massvacc(round(t),i)*pop(i+5*al); %dC/dt
        
        % Vaccine compartments
        dydt(i+6*al,1)=dydt(i+6*al,1) + veff*massvacc(round(t),i)*pop(i); % dV1/dt
        dydt(i+7*al,1)=dydt(i+7*al,1) + veff*massvacc(round(t),i)*(pop(i+2*al) + pop(i+3*al)); % dV2/dt 
        
        dydt(i+8*al,1)=dydt(i+8*al,1) + (1-veff)*massvacc(round(t),i)*pop(i); % dS1v/dt 
        dydt(i+9*al,1)=dydt(i+9*al,1) + massvacc(round(t),i)*pop(i+al); % dIv/dt
        dydt(i+10*al,1)=dydt(i+10*al,1) + (1-veff)*massvacc(round(t),i)*pop(i+2*al); % dRv/dt
        dydt(i+11*al,1)=dydt(i+11*al,1) + (1-veff)*massvacc(round(t),i)*pop(i+3*al); % dS2v/dt
        dydt(i+12*al,1)=dydt(i+12*al,1) + massvacc(round(t),i)*pop(i+4*al); % dI2v/dt
        dydt(i+13*al,1)=dydt(i+13*al,1) + massvacc(round(t),i)*pop(i+5*al); % dCv/dt
        
        % Number of vaccines distributed during initial campaigns        
        dydt(i+15*al,1) = massvacc(round(t),i)*(pop(i)+pop(i+al)+pop(i+2*al)+pop(i+3*al)+pop(i+4*al)+pop(i+5*al)); % dydt(i+15*al,1) +
    end

end
  
end

