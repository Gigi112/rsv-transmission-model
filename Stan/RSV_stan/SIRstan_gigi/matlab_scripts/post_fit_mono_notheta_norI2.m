% Post Stan fits
% This is the code we used for the publication

% cd('/Users/Marina/Dropbox/SIRstan_mono');
cd('C:\Users\marina_rdp\Documents\Dropbox\SIRstan_mono\')

% cd('/Users/Marina/Dropbox/SIRstan_mono')

output.kolkata = readtable('./kolkata/outjuly3.csv');
output.delhi = readtable('./delhi/outjuly3.csv');
output.kibera = readtable('./kibera/outjuly3.csv');
output.lwak = readtable('./lwak/outjuly3.csv');
output.dong_thap = readtable('./dong_thap/outjuly3.csv');

output.kolkata = output.kolkata(1:5:end,:);
output.delhi = output.delhi(1:5:end,:);
output.kibera = output.kibera(1:5:end,:);
output.lwak = output.lwak(1:5:end,:);
output.dong_thap = output.dong_thap(1:5:end,:);

sites = {'kolkata', 'delhi', 'kibera',  'lwak', 'dong_thap'}; % 'jakarta', 'dhaka', 'lwak',
Sites = {'Kolkata', 'Delhi', 'Nairobi',  'Lwak', 'Dong Thap'};
% for i=1:length(sites)
% beta.(sites{i}) = output.(sites{i}).beta;
% prop.(sites{i}) = output.(sites{i}).prop;
% mult1.(sites{i}) = exp(-output.(sites{i}).logm1);
% mult2.(sites{i}) = 1-exp(-output.(sites{i}).logm2);
% sensrate.(sites{i}) = output.(sites{i}).sensrate;
% r.(sites{i}) = output.(sites{i}).r;
% end

% Pull in fixed parameters

rng('shuffle')
% seed = RandStream.getGlobalStream.Seed;

yr_br = readtable('./data/birth_rates.csv'); % MULTIPLY TIMES POP?

data.kolkata = readtable('./data/kolkata.csv');
data.jakarta = readtable('./data/jakarta.csv');
data.kibera = readtable('./data/kibera.csv');
data.lwak = readtable('./data/lwak.csv');
data.dhaka = readtable('./data/dhaka.csv');
data.naheed = readtable('./data/naheed.csv');
data.delhi = readtable('./data/delhi.csv');
data.dong_thap = readtable('./data/dong_thap.csv');

% Constant parameters
params.delta = 1/4;
params.alpha = 0.01; % 0.01;
params.omega = 1/104; % -log(1-.25)/52;
% params.epsilon = 2.68*10^(-11);

% Deal with birthrates
wbr.kolkata = log(1+table2array(yr_br(1,2))*2)/52;
wbr.jakarta = log(1+table2array(yr_br(2,2)))/52;
wbr.lwak = log(1+table2array(yr_br(3,2)))/52;
wbr.kibera = log(1+table2array(yr_br(4,2)))/52;
wbr.dhaka = log(1+table2array(yr_br(5,2)))/52;
wbr.delhi = log(1+table2array(yr_br(7,2)))/52;
wbr.dong_thap = log(1+table2array(yr_br(8,2)))/52;

for i = 1:length(sites)
% study = sites{i};
agepar.(sites{i}).u = data.(sites{i}).aging;
% agepar.(sites{i}).mu_R0 = -log(1-data.(sites{i}).yr_mortality)/52;
% agepar.(sites{i}).mu(end) = agepar.(sites{i}).u(end-1)*data.(sites{i}).pyo(end-1)/data.(sites{i}).pyo(end);

agepar.(sites{i}).theta = data.(sites{i}).chronic;
agepar.(sites{i}).theta2 = zeros(length(agepar.(sites{i}).theta) ,1);% eval([study, '.chronic']);
agepar.(sites{i}).mub = [wbr.(sites{i}); zeros(length(agepar.(sites{i}).u)-1, 1)];

agepar.(sites{i}).participation = data.(sites{i}).draw;
agepar.(sites{i}).under5 = data.(sites{i}).under5; % under 5 or not?

agepar.(sites{i}).vol = data.(sites{i}).vol; 

cases.(sites{i}) = data.(sites{i}).cases;
population.(sites{i}) = data.(sites{i}).pyo;

agepar.(sites{i}).mu = ([wbr.(sites{i}); agepar.(sites{i}).u(1:(end-1), 1)].*[sum(data.(sites{i}).pyo); population.(sites{i})(1:(end-1),1)])./population.(sites{i})(:,1); 
agepar.(sites{i}).mu = agepar.(sites{i}).mu - agepar.(sites{i}).u;
end

% Sample from output.

cd('./matlab_scripts')

% iterations = 5000;
% draw=100;
randout = randsample(1:length(output.delhi.beta),draw);

% for i=1:length(sites)
% 
%         randnow = randout(j);
%         estpar.beta = output.(sites{i}).beta;
%         estpar.mult1 = exp(-output.(sites{i}).logm1);
%         estpar.mult2 = 1-exp(-output.(sites{i}).logm2);
%         estpar.prop = output.(sites{i}).prop;
%         estpar.r = 1; % r.(sites{i})(randi);
%         estpar.rC = output.(sites{i}).r;
%         params.epsilon = 0;
%         estpar.sensrate = output.(sites{i}).sensrate(randnow);
%         blood = 1-exp(-agepar.(sites{i}).vol*estpar.sensrate);
        
tic
for i=1:length(sites)
    popdist = population.(sites{i})./sum(population.(sites{i}));
    for j=1:draw
        randnow = randout(j);
        estpar.beta = output.(sites{i}).beta(randnow);
        estpar.mult1 = exp(-output.(sites{i}).logm1(randnow));
        estpar.mult2 = 1-exp(-output.(sites{i}).logm2(randnow));
        estpar.prop = output.(sites{i}).prop(randnow);
        estpar.r = 1; % r.(sites{i})(randi);
        estpar.rC = output.(sites{i}).r(randnow);
        params.epsilon = 0;

       % Simulation equation
       if i==4
        tspan=200*52+1;
       else
        tspan=100*52+1;
       end
        [postfit.(sites{i}).eqmcases, postfit.(sites{i}).eqmchronic, postfit.(sites{i}).eqmpop] = fn_forward_sim(params, agepar.(sites{i}), estpar, population.(sites{i}), tspan, sites{i}); 

        % Then adjust for the observation process.
        estpar.sensrate = output.(sites{i}).sensrate(randnow);
        blood = 1-exp(-agepar.(sites{i}).vol*estpar.sensrate);
        postfiteqmcases.(sites{i})(j,:) = postfit.(sites{i}).eqmcases;
        postfiteqmchronic.(sites{i})(j,:) = postfit.(sites{i}).eqmchronic;
        postfiteqmpop.(sites{i})(j,:) = postfit.(sites{i}).eqmpop;
        postfitrep.(sites{i})(j,:) = [(postfit.(sites{i}).eqmcases).*blood.*agepar.(sites{i}).participation*estpar.prop]';
        
        newbeta = [estpar.mult1*estpar.beta; estpar.mult2*estpar.beta; estpar.beta*ones(length(popdist)-2,1)];
        % R0 = beta/(mu+delta)*(1+r*theta*delta/mu)
        R0_samp(j,i) = newbeta'*popdist/(agepar.(sites{i}).mu'*popdist+params.delta)*(1+(estpar.rC*(agepar.(sites{i}).theta'*popdist)*params.delta)./(agepar.(sites{i}).mu'*popdist));
    end
end
toc

% More R0, without running the ODE...
% RUN R0 with mub (birth rates) not mu (death rates)
% Can't do this as a matrix? So I can have a full vector of R0?
% Should code it into stan.

% draw=1000;
% tic
% for i=1:length(sites)
%     popdist = population.(sites{i})./sum(population.(sites{i}));
%     for j=1:draw
%         randnow = randsample(1:length(output.(sites{i}).logm1),1);
%         estpar.beta = output.(sites{i}).beta(randnow);
%         estpar.mult1 = exp(-output.(sites{i}).logm1(randnow));
%         estpar.mult2 = 1-exp(-output.(sites{i}).logm2(randnow));
%         estpar.prop = output.(sites{i}).prop(randnow);
%         estpar.r = 1; % r.(sites{i})(randi);
%         estpar.rC = output.(sites{i}).r(randnow);
%         params.epsilon = 0;
% 
% 
%         newbeta = [estpar.mult1*estpar.beta; estpar.mult2*estpar.beta; estpar.beta*ones(length(popdist)-2,1)];
%         % R0 = beta/(mu+delta)*(1+r*theta*delta/mu)
%         R0_samp_more(j,i) = newbeta'*popdist/(agepar.(sites{i}).mu'*popdist+params.delta)*(1+(estpar.rC*(agepar.(sites{i}).theta'*popdist)*params.delta)./(agepar.(sites{i}).mu'*popdist));
%     end
% end
% toc

save('./samples_postfit_mono_notheta_norI2_jul20.mat', 'postfitrep', ...
    'postfiteqmcases', 'postfiteqmpop', 'R0_samp')
% , 'R0_samp_more'

% prctile(R0_samp_more, [2.5 50 97.5])';

% Graphs to look at the fits
% Clustered bar graphs with error bars

% figure('position', [0, 0, 1500, 300], 'color', 'w')
% [ha, pos] = tight_subplot(1, 5,[.025 .025],[.2 .2],[.05 .05]); 

figure('position', [0, 0, 1000, 600], 'color', 'w')
[ha, pos] = tight_subplot(2, 3,[.125 .05],[.1 .1],[.075 .075]); 
% gap, marg_h, marg_w

colormap([0.2, 0.2, 0.2; 0.8, 0.8, 0.8])
j=[1,2,4,5,3];

for i=[1, 2, 5, 3, 4]
    clear Y errY 
    Y = [cases.(sites{i})./population.(sites{i}), (median(postfitrep.(sites{i})./postfiteqmpop.(sites{i}),1)')]*1e5;
    errY(:,:,1) = [cases.(sites{i})./population.(sites{i}) - gaminv(0.025, cases.(sites{i}), 1./population.(sites{i})), [median(postfitrep.(sites{i})./postfiteqmpop.(sites{i}), 1) - prctile(postfitrep.(sites{i})./postfiteqmpop.(sites{i}), 2.5)]']*1e5; 
    errY(:,:,2) = [gaminv(0.975, max(cases.(sites{i}),1), 1./population.(sites{i})) - cases.(sites{i})./population.(sites{i}), [prctile(postfitrep.(sites{i})./postfiteqmpop.(sites{i}), 97.5) - median(postfitrep.(sites{i})./postfiteqmpop.(sites{i}), 1)]']*1e5;

    gaminv([0.025, 0.975], 1,1/(1*1000+1))

    axes(ha(j(i)))
    barwitherr(errY, Y)
    xlabel('Age Groups', 'FontSize',12)
    set(gca, 'XTickLabel', data.(sites{i}).age_group, 'XTickLabelRotation', 45, 'FontSize',12) % 'yscale', 'log'
    title(Sites{i}, 'FontSize',14)
    box(ha(i), 'on')
end

axes(ha(6))
barwitherr(zeros(2, 2, 2), zeros(2,2), 'Showbaseline', 'off', 'EdgeColor', 'w')
axis off
legend({'Observed Incidence', 'Model-Predicted Incidence'}, 'Location', [0.7, 0.25, 0.18, 0.1], 'FontSize', 12)

% [ax1,h1] = suplabel('Observed and model-predicted incidence', 't', [0.0100 0.1600 0.9800 0.7700])
% set(h1,'FontSize',16)
[ax2,h2] = suplabel('Cases per 100,000 people', 'y', [0.060 0.0600 0.9300 0.8800])
set(h2,'FontSize',14)

% MUST add a legend...

figure('position', [0, 0, 600, 700], 'color', 'w')
[ha, pos] = tight_subplot(3,2,[.05 .05],[.1 .1],[.1 .1]);

for i=1:length(sites)
    clear Y errY 
    Y = [median(postfiteqmchronic.(sites{i})./postfiteqmpop.(sites{i}),1)'];
    errY(:,:,1) = [median(postfiteqmchronic.(sites{i})./postfiteqmpop.(sites{i}), 1) - prctile(postfiteqmchronic.(sites{i})./postfiteqmpop.(sites{i}), 2.5)]'; 
    errY(:,:,2) = [prctile(postfiteqmchronic.(sites{i})./postfiteqmpop.(sites{i}), 97.5) - median(postfiteqmchronic.(sites{i})./postfiteqmpop.(sites{i}), 1)]';

    axes(ha(i))
    barwitherr(errY, Y)
    % set(gca, 'yscale', 'log') 
    box(ha(i), 'on')
end

figure
for i=1:length(sites)
    allchronics(:,i) = sum(postfiteqmchronic.(sites{i}),2)./sum(postfiteqmpop.(sites{i}),2);
end

chronics = median(allchronics)';
err_chronics(:,1) = prctile(allchronics, 50) - prctile(allchronics, 2.5); 
err_chronics(:,2) = prctile(allchronics, 97.5) - prctile(allchronics, 50);
barwitherr(err_chronics, chronics)


%% Vaccine simulations...

% Simulate the same age groups for all places this time...

agevac_label = {'0-<9m', '9m-<2y', '2y-<5y', '5-<10', '10-<15', ...
            '15-<20', '20-<25', '25-<30', '30-<40', '40-<50', '50+'}; 

% Vaccine effect parameters
% load('/Users/Marina/Google Drive/Dissertation/Vax Eff Estimate/vax_samp_July2016.mat')
load('C:/Users/Marina Antillon/Google Drive/Dissertation/Vax Eff Estimate/vax_samp_July2016.mat')
% column #1 is probability of 
% Load data_vac, parameters for each place for a 'standard 
% population as above...'

data.kolkata = readtable('../data/kolkata_vax.csv');
data.kibera = readtable('../data/kibera_vax.csv');
data.lwak = readtable('../data/lwak_vax.csv');
data.delhi = readtable('../data/delhi_vax.csv');
data.dong_thap = readtable('../data/dong_thap_vax.csv');

clear agepar
for i = 1:length(sites)
population.(sites{i}) = data.(sites{i}).pyo;
agepar.(sites{i}).u = data.(sites{i}).aging;

agepar.(sites{i}).theta = [ones(6,1)*0.003; 0.021; 0.021; 0.044; 0.088; 0.101];
agepar.(sites{i}).theta2 = zeros(length(agepar.(sites{i}).theta) ,1);% eval([study, '.chronic']);
agepar.(sites{i}).mub = [wbr.(sites{i}); zeros(length(agepar.(sites{i}).theta)-1, 1)];
agepar.(sites{i}).under5 = data.(sites{i}).under5; % under 5 or not?

agepar.(sites{i}).mu = ([wbr.(sites{i}); agepar.(sites{i}).u(1:(end-1), 1)].*[sum(data.(sites{i}).pyo); population.(sites{i})(1:(end-1),1)])./population.(sites{i})(:,1); 
agepar.(sites{i}).mu = agepar.(sites{i}).mu - agepar.(sites{i}).u;

% agepar.(sites{i}).mu_real = -log(1-data.(sites{i}).yr_mortality)/52;
agepar.(sites{i}).mu_real = agepar.(sites{i}).mu;
end

% May or may not be necessary for this portion of the analysis.
% yrlived = [4.5/12, (9+24)/2, 3.5, 7.5, 12.5, 17.5, 22.5, 27.5, 35, 45, nan]; 
% For the last category, consider the life expectancy in each place...
% aging_vac = [1/(9/12), 1/((24-9)/12), 1/3, 1/5*ones(1,5), 1/10*ones(1,2), 0]./52.25;

vcov1 = 0.8;
vcov2 = 0.7;
        
int_plans = [1, 1, 0, 0; ... % No intervention
    2, 2, vcov1, 0; ... % Routine vaccination (the 3 in the second position doesn't matter)
    2, 3, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 5 years of age
    2, 5, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 15 years of age
    2, 7, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 25 years of age
    2, 11, vcov1, vcov2]; % ; ...  % Routine vaccination and catch up campaign for everyone
    % 4, 4, vcov1, 0; ... % Routine vaccination at 5 years of age
    % 4, 5, vcov1, vcov2]; % Routine vaccination at 5 years of age and catch up campaign for people up to 15 years of age

draw=1000;
randout = randsample(1:length(output.delhi.beta),draw); 
% if draw=length(output.delhi.beta) then it is just using all of the 
% draws from stan. (this is sampling without replacement)
randout2 = randsample(1:length(vax_samp.nu),draw);

% the list of simulationsnames must be in the same order as the details of
% the campaigns in the int_plans matrix above.
simulationnames = {'Base_site', 'Routine_site',  'RoutineC05_site', ...
                   'RoutineC15_site', 'RoutineC25_site', 'RoutineC99_site'...
                   'Routine05_site', 'Routine05C15_site'};

% Make population 100,000. Simulation #1 runs with mu_equilibrium,
% simulation #2 (w vaccination) runs for 30 years with actual mu.

% clear 'Base' 'Routine' 'RoutineC05' 'RoutineC15'
% clear 'RoutineC25' 'RoutineC99' 'ParamsVax' 'ParamsDis'

% myCluster = parcluster('local')
% delete(myCluster.Jobs)
delete(gcp) % Terminating any existing sessions
% isempty(gcp('nocreate'))
% parpool('local', 5) % connecting to 4 workers

parpool(length(sites))

tic
clear Base Routine RoutineC05 RoutineC15 RoutineC25 RoutineC99 Routine05 Routine05C15 ParamsDis ParamsVax
parfor i=1:length(sites)

[Base_site, Routine_site, RoutineC05_site, RoutineC15_site, RoutineC25_site, RoutineC99_site, ParamsDis_site, ParamsVax_site] = fn_vax_sim(output, vax_samp, agepar, params, randout, randout2, population, int_plans, simulationnames, sites, draw, i);
    %  Routine05_site, Routine05C15_site, 
    Base(:,i) = Base_site;
    Routine(:,i) = Routine_site;
    RoutineC05(:,i) = RoutineC05_site;
    RoutineC15(:,i) = RoutineC15_site;
    RoutineC25(:,i) = RoutineC25_site;
    RoutineC99(:,i) = RoutineC99_site;
    
    % Routine05(:,i) = Routine05_site;
    % Routine05C15(:,i) = Routine05C15_site;

    ParamsVax(:,i) = ParamsVax_site;
    ParamsDis(:,i) = ParamsDis_site;
    
end
toc

save('vax_sim_Nov21.mat', 'Base', 'Routine',  'RoutineC05', 'RoutineC15', 'RoutineC25', 'RoutineC99', 'ParamsVax', 'ParamsDis')
% 'Routine05', 'Routine05C15',

%% Vaccine simulations WITH water

% Simulate the same age groups for all places this time...

agevac_label = {'0-<9m', '9m-<2y', '2y-<5y', '5-<10', '10-<15', ...
            '15-<20', '20-<25', '25-<30', '30-<40', '40-<50', '50+'}; 

% Vaccine effect parameters
% load('/Users/Marina/Google Drive/Dissertation/Vax Eff Estimate/vax_samp_July2016.mat')
load('C:/Users/Marina Antillon/Google Drive/Dissertation/Vax Eff Estimate/vax_samp_July2016.mat')
% column #1 is probability of 
% Load data_vac, parameters for each place for a 'standard 
% population as above...'

data.kolkata = readtable('../data/kolkata_vax.csv');
data.kibera = readtable('../data/kibera_vax.csv');
data.lwak = readtable('../data/lwak_vax.csv');
data.delhi = readtable('../data/delhi_vax.csv');
data.dong_thap = readtable('../data/dong_thap_vax.csv');

clear agepar
for i = 1:length(sites)
population.(sites{i}) = data.(sites{i}).pyo;
agepar.(sites{i}).u = data.(sites{i}).aging;

agepar.(sites{i}).theta = [ones(6,1)*0.003; 0.021; 0.021; 0.044; 0.088; 0.101];
agepar.(sites{i}).theta2 = zeros(length(agepar.(sites{i}).theta) ,1);% eval([study, '.chronic']);
agepar.(sites{i}).mub = [wbr.(sites{i}); zeros(length(agepar.(sites{i}).theta)-1, 1)];
agepar.(sites{i}).under5 = data.(sites{i}).under5; % under 5 or not?

agepar.(sites{i}).mu = ([wbr.(sites{i}); agepar.(sites{i}).u(1:(end-1), 1)].*[sum(data.(sites{i}).pyo); population.(sites{i})(1:(end-1),1)])./population.(sites{i})(:,1); 
agepar.(sites{i}).mu = agepar.(sites{i}).mu - agepar.(sites{i}).u;

% agepar.(sites{i}).mu_real = -log(1-data.(sites{i}).yr_mortality)/52;
agepar.(sites{i}).mu_real = agepar.(sites{i}).mu;
end

% May or may not be necessary for this portion of the analysis.
% yrlived = [4.5/12, (9+24)/2, 3.5, 7.5, 12.5, 17.5, 22.5, 27.5, 35, 45, nan]; 
% For the last category, consider the life expectancy in each place...
% aging_vac = [1/(9/12), 1/((24-9)/12), 1/3, 1/5*ones(1,5), 1/10*ones(1,2), 0]./52.25;

vcov1 = 0.8;
vcov2 = 0.7;
        
int_plans = [1, 1, 0, 0; ... % No intervention
    2, 2, vcov1, 0; ... % Routine vaccination (the 3 in the second position doesn't matter)
    2, 3, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 5 years of age
    2, 5, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 15 years of age
    2, 7, vcov1, vcov2; ... % Routine vaccination and catch up campaign for people up to 25 years of age
    2, 11, vcov1, vcov2]; % ; ...  % Routine vaccination and catch up campaign for everyone
    % 4, 4, vcov1, 0; ... % Routine vaccination at 5 years of age
    % 4, 5, vcov1, vcov2]; % Routine vaccination at 5 years of age and catch up campaign for people up to 15 years of age

draw=500;
randout = randsample(1:length(output.delhi.beta),draw); 
% if draw=length(output.delhi.beta) then it is just using all of the 
% draws from stan. (this is sampling without replacement)
randout2 = randsample(1:length(vax_samp.nu),draw);


% the list of simulationsnames must be in the same order as the details of
% the campaigns in the int_plans matrix above.
simulationnames = {'Base_site', 'Routine_site',  'RoutineC05_site', ...
                   'RoutineC15_site', 'RoutineC25_site', 'RoutineC99_site'};
% 'Routine05_site', 'Routine05C15_site'
% Make population 100,000. Simulation #1 runs with mu_equilibrium,
% simulation #2 (w vaccination) runs for 30 years with actual mu.

% clear 'Base' 'Routine' 'RoutineC05' 'RoutineC15'
% clear 'RoutineC25' 'RoutineC99' 'ParamsVax' 'ParamsDis'

% myCluster = parcluster('local')
% delete(myCluster.Jobs)
delete(gcp) % Terminating any existing sessions
% isempty(gcp('nocreate'))
% parpool('local', 5) % connecting to 4 workers

parpool(length(sites))
% Routine05 Routine05C15
tic
clear Base Routine RoutineC05 RoutineC15 RoutineC25 RoutineC99  ParamsDis ParamsVax
parfor i=1:length(sites)
tic
[Base_site, Routine_site, RoutineC05_site, RoutineC15_site, RoutineC25_site, RoutineC99_site, ParamsDis_site, ParamsVax_site] = fn_vax_sim_ext_h2o(output, vax_samp, agepar, params, randout, randout2, population, int_plans, simulationnames, sites, draw, i);
toc   
%  Routine05_site, Routine05C15_site, 
    Base(:,i) = Base_site;
    Routine(:,i) = Routine_site;
    RoutineC05(:,i) = RoutineC05_site;
    RoutineC15(:,i) = RoutineC15_site;
    RoutineC25(:,i) = RoutineC25_site;
    RoutineC99(:,i) = RoutineC99_site;
    
    % Routine05(:,i) = Routine05_site;
    % Routine05C15(:,i) = Routine05C15_site;

    ParamsVax(:,i) = ParamsVax_site;
    ParamsDis(:,i) = ParamsDis_site;
    
end
toc

save('vax_sim_Apr25_h2o.mat', 'Base', 'Routine',  'RoutineC05', 'RoutineC15', 'RoutineC25', 'RoutineC99', 'ParamsVax', 'ParamsDis')
% 'Routine05', 'Routine05C15',
%% This was last done with 'vax_sim_Nov21'
plot(Base(1,1).cumI1(2:31,:)-Base(1,1).cumI1(1:30,:))
bar([Base(1,i).cumI1(2:31,:)-Base(1,i).cumI1(1:30,:)]'./repmat(Base(1,i).cumI1(2,:)',1,30))
plot(Base(1,1).pop(2:31,:))

PDE = nan(1000, 5);
PDE_C05 = nan(1000, 5);
PDE_C15 = nan(1000, 5);
PDE_C25 = nan(1000, 5);
PDE_C99 = nan(1000, 5);

for i = 1:1000
    for j = 1:5
        Dnovacc = Base(i,j).cumI1(2:11,:)-Base(i,j).cumI1(1:10,:);
        
        Vprop = (Routine(i,j).V1(2:11,:)+Routine(i,j).V2(2:11,:))./Routine(i,j).pop(2:11,:);
        Ddynamic = sum(sum((Routine(i,j).cumI1(2:11,:)-Routine(i,j).cumI1(1:10,:))));
        PDE(i,j) = 1-(sum(sum(Dnovacc.*(1-Vprop)))-Ddynamic)/(sum(sum(Dnovacc))-Ddynamic);
    
        Vprop = (RoutineC05(i,j).V1(2:11,:)+RoutineC05(i,j).V2(2:11,:))./RoutineC05(i,j).pop(2:11,:);
        Ddynamic = sum(sum((RoutineC05(i,j).cumI1(2:11,:)-RoutineC05(i,j).cumI1(1:10,:))));
        PDE_C05(i,j) = 1-(sum(sum(Dnovacc.*(1-Vprop)))-Ddynamic)/(sum(sum(Dnovacc))-Ddynamic);
        
        Vprop = (RoutineC15(i,j).V1(2:11,:)+RoutineC15(i,j).V2(2:11,:))./RoutineC15(i,j).pop(2:11,:);
        Ddynamic = sum(sum((RoutineC15(i,j).cumI1(2:11,:)-RoutineC15(i,j).cumI1(1:10,:))));
        PDE_C15(i,j) = 1-(sum(sum(Dnovacc.*(1-Vprop)))-Ddynamic)/(sum(sum(Dnovacc))-Ddynamic);

        Vprop = (RoutineC25(i,j).V1(2:11,:)+RoutineC25(i,j).V2(2:11,:))./RoutineC25(i,j).pop(2:11,:);
        Ddynamic = sum(sum((RoutineC25(i,j).cumI1(2:11,:)-RoutineC25(i,j).cumI1(1:10,:))));
        PDE_C25(i,j) = 1-(sum(sum(Dnovacc.*(1-Vprop)))-Ddynamic)/(sum(sum(Dnovacc))-Ddynamic);

        Vprop = (RoutineC99(i,j).V1(2:11,:)+RoutineC99(i,j).V2(2:11,:))./RoutineC99(i,j).pop(2:11,:);
        Ddynamic = sum(sum((RoutineC99(i,j).cumI1(2:11,:)-RoutineC99(i,j).cumI1(1:10,:))));
        PDE_C99(i,j) = 1-(sum(sum(Dnovacc.*(1-Vprop)))-Ddynamic)/(sum(sum(Dnovacc))-Ddynamic);
    end
end

% Remember that Vietnam is last here, and it is in the middle elsewhere
pdestruct = struct([]);
pdestruct{1,1} = prctile(PDE, [0.5 0.025 0.975]);
pdestruct{2,1} = prctile(PDE_C05, [0.5 0.025 0.975]);
pdestruct{3,1} = prctile(PDE_C15, [0.5 0.025 0.975]);
pdestruct{4,1} = prctile(PDE_C25, [0.5 0.025 0.975]);
pdestruct{5,1} = prctile(PDE_C99, [0.5 0.025 0.975]);

save('pdestruct.mat', 'pdestruct')

%% OLD SCRAPS

% It seems we need to run Lwak for 100-150 years to really see it reach
% equilibrium...

% to check population changes:
% bar([vax04(1,1).eqmpop./repmat(sum(vax04(1,1).eqmpop,2), 1, 11)], 'stacked')


% Dimensions: 5 places (struct/list level 1), 4 interventions (struct/list level 2), 
% 11 age groups, 31 years, 1000 iterations, P measures...

% Save a file with age groups, but sum over all age groups before sending
% the file to R.

% How to pull out vaccine effects... ODE gives incidence, multiply by
% original population.
% Pull out chronic carrier population...
% Pull out vaccination rates
% Aggregate things by year, but keep year-to-year effects.
% Remember we are taking into account cases even if the blood 
    % is not drawn or the culture does not turn out positive.

% Interventions (all ViCV)
% Assuming how much coverage.
% None, R 9m, R 9m + C to 5y, R 9m + C to 15y, 9m + C to 25y, 9m + M
% Do second option here instead of 9m, 2y? Because the characteristics of
% the vacc?
% R 6y, R 6y + C to 15y, R 6y + C to 25y, R 6y + M 

% Separately: hospitalization rates, deaths, etc. all in R, 

% Condense the output.

simulationnames = {'Base', 'Routine',  'RoutineC05', 'RoutineC15', 'RoutineC25', 'RoutineC99'};

% - median baseline incidence (with reporting rate, but not BCS)
sum(Base(i, k).cumI1,2)

% - median direct effect (incidence among cumIv only).
sum(Routine(i, k).cumI1v,2) % How many people were vaccinated? 
% Going to need to recalculate it...

% - median, LL, UL overall effect (incidence among cumI1v+cumI1)
OE.(sites{i})(iteration, 1:31) = sum(Routine(i, k).cumI1,2)' + sum(Routine(i, k).cumI1v,2)';

% Grouped stacked bars
% Group is site, bar is intervention, color is age group: <5, 5-15, 15-25, 25+
% Median, LL, UL overall effect  

horizon = 10;

for intervention = 1:6
    for k = 1:size(Base, 2)
        eval(['tmp = ' simulationnames{intervention} ';']);
        for i = 1:size(Base, 1)
            OE.(sites{k}).OEage1(i, intervention) = sum(sum(tmp(i, k).cumI1(1:3,1:horizon+1),2)' + sum(tmp(i, k).cumI1v(1:3,1:horizon+1),2)');
            OE.(sites{k}).OEage2(i, intervention) = sum(sum(tmp(i, k).cumI1(4:5,1:horizon+1),2)' + sum(tmp(i, k).cumI1v(4:5,1:horizon+1),2)');
            OE.(sites{k}).OEage3(i, intervention) = sum(sum(tmp(i, k).cumI1(6:7,1:horizon+1),2)' + sum(tmp(i, k).cumI1v(6:7,1:horizon+1),2)');
            OE.(sites{k}).OEage4(i, intervention) = sum(sum(tmp(i, k).cumI1(9:11,1:horizon+1),2)' + sum(tmp(i, k).cumI1v(9:11,1:horizon+1),2)');
        end
    end
end

% dimensions: , 
% chronics(,) = median(OE.(sites{k}).OEage1(i, intervention))', ;

err_chronics(:,1) = prctile(allchronics, 50) - prctile(allchronics, 2.5); 
err_chronics(:,2) = prctile(allchronics, 97.5) - prctile(allchronics, 50);
barwitherr(err_chronics, chronics)

% 10 year horizon... (maybe in supplement give a 20-30 year horizon).

% Graph the vaccine interventions...
figure
% loop over the interventions
for j=1:4
subplot(2,4,j)
hold on

fill([1/52:1/52:(tspan-t0)/52 (tspan-t0)/52:-1/52:1/52 1/52]',52.18*100000*[overallimpact_dr1_lCI(:,j); flipud(overallimpact_dr1_uCI(:,j)); overallimpact_dr1_lCI(1,j)],[1 .8 .8],'EdgeColor',[1 .8 .8])
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact_dr1_lCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*overallimpact_dr1_uCI(:,j),'--r','LineWidth',1)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*(sum(vaccimpact_dr1(1,j).Dnovacc,2)-vaccimpact_dr1(1,j).popdirecteff)./sum(vaccimpact_dr1(1,j).pop,2),'g','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).D,2)./sum(vaccimpact_dr1(1,j).pop,2),'r','LineWidth',2)
plot(1/52:1/52:(tspan-t0)/52,52.18*100000*sum(vaccimpact_dr1(1,j).Dnovacc,2)./sum(vaccimpact_dr1(1,j).pop,2),'--b','LineWidth',2)
%ylim([0 250]) %ulim(j)])
if j==1
ylabel({'Incidence rate'; '(per 100,000 per year)'},'FontSize',12) %,'FontSize',14)
elseif j==2
legend('No vaccination','Direct effect of vaccination','Predicted overall effect','Location','SO','Orientation','Horizontal')
%ylim([0 800])
end
xlabel('Year') %,'FontSize',14)
title(inciddata2(j).setting,'FontWeight','Bold','FontSize',12)

subplot(2,4,4+j)
hold on
bar(vaccimpact_dr1(1,j).cumincidnovacc(11,:)-vaccimpact_dr1(1,j).cumincid(11,:))
plot([1:3; 1:3],[prctile(casesaverted_dr1_10y(:,:,j),2.5); prctile(casesaverted_dr1_10y(:,:,j),97.5)],'k')
set(gca,'XTick',1:3,'XTickLabel',{'<5y','5-20y','20+y'})
xlabel('Age group')
if j==1
ylabel({'Cases averted (per 100,000)';'over 10 years'},'FontSize',12)
end
end