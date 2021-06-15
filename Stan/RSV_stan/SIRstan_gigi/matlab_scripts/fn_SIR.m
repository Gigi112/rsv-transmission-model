function dydt = fn_SIR(t, pop, beta, mub, u, mu, omega, r, rC, alpha, delta, theta, theta2, epsilon, al, population)
%Differential equations for typhoid model.

% make two more parameters, and then we are ready for ode's
lambda = beta*sum(pop((al+1):(2*al))+r*pop((4*al+1):(5*al))+rC*pop((5*al+1):(6*al)))/sum(population);
births = mub*sum(population);

for i=1:al
    dydt(i,1) = births(i) + epsilon*pop(i+3*al) - lambda(i)*pop(i) - (u(i)+mu(i))*pop(i); %dS1/dt
    dydt(i+al,1) = lambda(i)*pop(i) - delta*pop(i+al) - (u(i)+mu(i))*pop(i+al); %dI1/dt
    dydt(i+2*al,1) = delta*(1-alpha-theta(i))*pop(i+al) + delta*(1-theta2(i))*pop(i+4*al) - omega*pop(i+2*al) - (u(i)+mu(i))*pop(i+2*al); %dR/dt
    dydt(i+3*al,1) = omega*pop(i+2*al) - epsilon*pop(i+3*al) - lambda(i)*pop(i+3*al) - (u(i)+mu(i))*pop(i+3*al); %dS2/dt 
    dydt(i+4*al,1) = lambda(i)*pop(i+3*al) - delta*pop(i+4*al) - (u(i)+mu(i))*pop(i+4*al); %dI2/dt
    dydt(i+5*al,1) = delta*(theta(i)*pop(i+al) + theta2(i)*pop(i+4*al)) - (u(i)+mu(i))*pop(i+5*al); %dC/dt
    
    dydt(i+6*al,1) = lambda(i)*pop(i); % Cumulative incidence of symptomatic infection, no one ages out of here...
    
    if i>1 %aging
        dydt(i,1)=dydt(i,1) + u(i-1)*pop(i-1);
        dydt(i+al,1)=dydt(i+al,1) + u(i-1)*pop(i-1+al);
        dydt(i+2*al,1)=dydt(i+2*al,1) + u(i-1)*pop(i-1+2*al);
        dydt(i+3*al,1)=dydt(i+3*al,1) + u(i-1)*pop(i-1+3*al);
        dydt(i+4*al,1)=dydt(i+4*al,1) + u(i-1)*pop(i-1+4*al); 
        dydt(i+5*al,1)=dydt(i+5*al,1) + u(i-1)*pop(i-1+5*al);
    end
end

end
