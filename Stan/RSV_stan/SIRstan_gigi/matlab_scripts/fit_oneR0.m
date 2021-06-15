% Bring in all data + parameters
clear
cd('/Users/Marina/Dropbox/SIRstan_mono')
% cd('C:\Users\marina_rdp\Documents\Dropbox\SIRstan_mono')

% Don't give me the same 'random' draws in every session.
rng('shuffle')
seed = RandStream.getGlobalStream.Seed;

yr_br = readtable('./data/birth_rates.csv');

kolkata = readtable('./data/kolkata.csv');
jakarta = readtable('./data/jakarta.csv');
kibera = readtable('./data/kibera.csv');
lwak = readtable('./data/lwak.csv');
dhaka = readtable('./data/dhaka.csv');
naheed = readtable('./data/naheed.csv');
delhi = readtable('./data/delhi.csv');
dong_thap = readtable('./data/dong_thap.csv');

% Deal with birthrates
wbr.kolkata = log(1+table2array(yr_br(1,2))*2)/52;
wbr.jakarta = log(1+table2array(yr_br(2,2)))/52;
wbr.lwak = log(1+table2array(yr_br(3,2)))/52;
wbr.kibera = log(1+table2array(yr_br(4,2)))/52;
wbr.dhaka = log(1+table2array(yr_br(5,2)))/52;
wbr.delhi = log(1+table2array(yr_br(7,2)))/52;
wbr.dong_thap = log(1+table2array(yr_br(8,2)))/52;

sites = {'kolkata', 'delhi', 'kibera', 'lwak', 'dong_thap'}; %'dhaka', 

% 'jakarta', 

for i = 1
study = sites{i};
site = sites{i};

clear data 

% Constant parameters
data.delta = 1/4;
data.alpha = 0.01;
data.omega = 1/104;
% data.epsilon = 2.68*10^(-11);

data.agegroups = length(eval([study, '.cases']));
data.u(:,1) = eval([study, '.aging']);
% LAST ONE SHOULD BE MUCH BIGGER. NEW_MU: MU*THE SIZE OF AGING IN/MU*SIZE
% OF SECOND TO OLDEST/SIZE OF OLDEST. WILL AFFECT CHRONIC CARRIAGE...
data.theta(:,1) = eval([study, '.chronic']);
data.mub(:,1) = wbr.(study);

data.participation = eval([study, '.draw']);
% data.under5(:,i) = [eval([study, '.under5']); zeros(9-data.agegroups(i), 1)]; % under 5 or not?

data.vol(:,1) = eval([study, '.vol']); 

data.cases(:,1) = round(eval([study, '.cases']));
data.y0(:,1) = round(eval([study, '.pyo']));

% old: data.mu(:,1) = -log(1-eval([study, '.yr_mortality']))/52;
data.mu(:,1) = ([data.mub; data.u(1:(end-1), 1)].*[sum(data.y0); data.y0(1:(end-1),1)])./data.y0; 
data.mu(:,1) = data.mu(:,1) - data.u(:,1);

data.casesdhaka2(1) = round(naheed.cases(1)/naheed.pyo(1)*sum(dhaka.pyo(1:2)));
data.casesdhaka2(2) = round(naheed.cases(2)/naheed.pyo(2)*sum(dhaka.pyo(3:4)));

if strcmp(sites{i}, 'dhaka')
    data.dhaka=1;
else
    data.dhaka=0;
end

% tspan = [50; 50; 50; 100; 50];
% iterations = 5000;

data.population = sum(data.y0)';
data.y0 = [data.y0; zeros(data.agegroups*7, 1)];
data.y0(data.agegroups+3) = 1;
data.y0(data.agegroups*6+3) = 1;
data.y0(data.agegroups*7+3) = 1./data.population;
data.y0index = reshape([1:data.agegroups*6], data.agegroups, 6);
data.mub = data.mub.*data.population';

odetime = 10*[50; 50; 50; 100; 50];
data.ts = odetime(i)*52;
data.tspan = 1:data.ts;

data.thetaindicator = 0;
data.rindicator=0;

init.r = 0.3; 

% clear data.under5
allbetas = [0.17; 0.5; 0.25; 0.15; 0.2]; 
% init.logbeta = log(allbetas(i));
data.y0(data.agegroups*7+3) = allbetas(i)./data.population;

init.R0 = 2;
init.logm1 = -log(0.12);
init.logm2 = -log(1-0.3);
logprop = -log(1-[0.2, 0.9, 0.3, 0.9, 0.3]);
% logprop(5:6) = -log([0.9, 0.9]);
init.logprop = logprop(i);
init.sensrate = 0.1747;

sirmodel = fileread('./SIR_typhi_one_sampleR0.stan');

clock

% mkdir(sites{i})
cd(sites{i})

fit = stan('model_code', sirmodel, 'model_name', 'SIR_typhi_notheta_norI2_apr16_2017', ...
            'data', data, 'init', init, 'warmup', 10, 'iter', 10, ...
            'file_overwrite', true, 'chains', 1, 'verbose', true);
		
cd('..')
		
end

out.beta_samp = fit.extract('permuted',false).beta;
out.prop_samp = fit.extract('permuted',false).logprop;
out.mult1_samp = fit.extract('permuted',false).logm1;
out.mult2_samp = fit.extract('permuted',false).logm2;

out.r_samp = fit.extract('permuted',false).r;
out.sens_samp = fit.extract('permuted',false).sensrate;


% fit2 = stan('fit', fit, 'warmup', 100, 'iter', 100, 'data', data, 'init', init);
% 
% save('nodr_feb28_stan.mat', 'fit', 'out')
% save('nodr_feb28_stan_parameters.mat', 'out')

% 
./SIR_typhi_nor_one sample num_warmup=100 num_samples=100 adapt delta=.8 init_buffer=50 term_buffer=45 window=5 algorithm=hmc engine=nuts max_depth=4 stepsize=1 stepsize_jitter=0 init=temp.init.R data file=temp.data.R output1 refresh=5

timenow = clock;
    sprintf('%i iterations performed at %i %i %i', j, timenow(4), timenow(5), timenow(6))
    end

tspan=500*52;
[postfit.(sites{i}).eqmcases, postfit.(sites{i}).eqmchronic, postfit.(sites{i}).eqmpop] = fn_forward_sim(params, agepar.(sites{i}), init, population.(sites{i}), tspan); 
