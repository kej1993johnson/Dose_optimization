% This script loads in the traj structure which contains the dose response
% for all conditions. It creates the trajsum structure which will contain a
% mean vector and U(t) vector for each different concentration tested. This
% will be used to calibrate to the model.
 close all; clear all; clc
 
 %% Load in data structure 
S = load('../out/trajfit.mat');
traj= S.traj;

%% Separate by N0
uniqdose = [];
doselist = [];
for i = 1:length(traj)
    if ~isempty(traj(i).dosenum)
    if traj(i).dosenum==1 || traj(i).dosenum == 0
    doselist =vertcat(doselist, traj(i).dose);
    end
    end
end
uniqdose= unique(doselist);
%% Make a new structure which combines each dox concentration
    % percapita growth rate 
    % variance of per capita growth rate
 
 % find groups by N0

 
 for i = 1:length(uniqdose)
    trajsum(i).Cdox = [];
    trajsum(i).Nmat = [];
    trajsum(i).nreps = 0;
    trajsum(i).tmat = [];
 end

 for i = 1:length(uniqdose) % number of unique seed numbers
    for j = 1:length(traj)
        date = {'8-16-18'};
        if contains(traj(j).date, date)  % only want data from this run
                if traj(j).dose == uniqdose(i)
                    trajsum(i).nreps = trajsum(i).nreps +1;
                    trajsum(i).Cdox = traj(j).dose;
                    trajsum(i).color = traj(j).color;
                    trajsum(i).tmat = horzcat(trajsum(i).tmat,traj(j).time);
                    trajsum(i).Nmat = horzcat(trajsum(i).Nmat, traj(j).rawN);
                    trajsum(i).tdose = traj(j).tdose;
                end
        end
    end
 end
 
 %% Plot the raw data
 for i = 1:length(trajsum)
     for j = 1: trajsum(i).nreps
         plot(trajsum(i).tmat(:,j), trajsum(i).Nmat(:,j), 'color', trajsum(i).color)
         hold on
     end
 end
 xlabel('time (hours)')
 ylabel('N(t)')
 title('N(t) for different single pulse treatments')
 %% Again, clean data for fitting...
 Nfin = 5e4;
 for i = 1:length(trajsum)
    N = trajsum(i).Nmat;
    t = trajsum(i).tmat;
    i0 = find(t(:,1)>trajsum(i).tdose,1,'first'); % I arbitrarily search for a maximum in the first 200 hours
    iend = find(N(:,1)>=Nfin,1, 'first');
    if ~isempty(iend)
    tfit = t(i0:iend,:)-t(i0, :); 
    Nfit = N(i0:iend, :);
    end
    if isempty(iend)
        tfit = t(i0:end, :)-t(i0, :); 
        Nfit = N(i0:end, :);
    end
    trajsum(i).tfit =tfit;
    trajsum(i).Nfit =Nfit;
end
  %% Plot the cleaned data
  figure;
 for i = 1:length(trajsum)
     for j = 1: trajsum(i).nreps
         plot(trajsum(i).tfit(:,j), trajsum(i).Nfit(:,j), 'color', trajsum(i).color)
         hold on
     end
 end
 xlabel('time (hours)')
 ylabel('N(t)')
 title('N(t) for different single pulse treatments')

%% Now find mean and standard deviation vectors

for i = 1:length(trajsum)
    trajsum(i).Nmean = mean(trajsum(i).Nfit,2);
    trajsum(i).tvec = trajsum(i).tfit(:,1);
    trajsum(i).Nstd = std(trajsum(i).Nfit,0,2);
end
% Plot the average data 
  figure;
 for i = 1:6%length(trajsum)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         plot(trajsum(i).tvec, trajsum(i).Nmean + trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - trajsum(i).Nstd, 'color', trajsum(i).color)

 end
 xlabel('time (hours)')
 ylabel('N(t)')
 title('N(t) for different single pulse treatments')
 %% Test forward model
% Set some arbitrary things as inputs

kdrug = 0.025;
k = 0.5;
dt = 1;
% input time vectors for each different dose response

for i = 1:length(trajsum)
    ttest = [];
    ttest = 1:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest)); 
end
%%
% set parameters of forward model
carcap = 2e4;
S0=2e3; % set initial conditions (assume all cells sensitive to start)
R0 = 0; 
Cdox = 75; % nm doxorubicin
rs = 0.015;
rr = 0.005;
ds = 0.001;
dr = 0;
alpha = 0.0001;

pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];


Nsr = [];
tsum = [];
tdrug = 1; % since we are monitoring first treatment only

for i = 1:6
    tvec = round(trajsum(i).tvec(2:end),0);
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    
    Nsr = vertcat(Nsr, Nsri);
end

figure;
plot(tsum, Nsr(:,1), 'b*')
xlabel('time (hours)')
ylabel('Model predicted response to pulse treat')
title('Test ability to generate model trajectories')
%% Bayesian fit for gs & K using untreated control

% The simplified  model looks like this
%
% $$ N(t) = N0*K*(1/(N0 + (K-N0)exp(-gst));
% 
% where:
%
% * $N_0$ is the initial cell number
% * gs is the sensitive cell growth rate 
% * K is the carrying capacity


% Define transforms 
% single exponential
pfxform1 = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform1 = @(phat)[1 1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output



sigma = trajsum(1).Nstd(2:end);

ydataf = trajsum(1).Nmean;
N0 = ydata(1);
ydata = ydataf(2:end);
ytimef = trajsum(1).tvec;
ytime = ytimef(2:end);

% Set up forward models, fit all three nested versions of model
modelfun1 = @(p)simmodel1(p, ytime, N0); % single exponential model with death  


    % INITIAL GUESSES BASED ON DATA

    gguess = (yfxform(ydata(end))-yfxform(ydata(end-5)))/(ytime(end)-ytime(end-5)); 
    % alter initial guesses to prevent NaNs and zeros
    Kguess = ydata(end);
    if isnan(gguess)
        gguess = 1e-5;
    end
    if isinf(gguess)
        gguess = 0.8;
    end
    if gguess <= 0
        gguess = 1e-5;
    end
    
    % Initial guess matrices
    theta1 = [gguess, Kguess]; % and k

    
    % Write log likelihood function based on assumption of normally
    % distributed sampling error
    
    % Goal: maximize the probability of the data given the model. Normpdf
    % will output a probability of the data (x- 1st argument), given the
    % mean(expectation, here the model), and the variance, at each time point. 
    % take the log and minimize the NLL to maximize likeihood of data given
    % the model
    
    % single exponential with carrying capacity
    loglikelihood1 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfun1(pbxform1(phat))), sigma)));
    % single exponential with growth
   
    % Write objective functions for each model
    objfun1 = @(phat)-loglikelihood1(phat); 
    phatbest1 = fminsearch(objfun1, pfxform1(theta1));
    
    pi = pbxform1(phatbest1);
    gs = pi(1);
    carcap = pi(2); 
    singexpmodel = simmodel1(pbxform1(phatbest1), ytime, N0);
    
    figure;
    plot(ytime, ydata, '*', 'LineWidth', 3)
    hold on
    plot(ytime, singexpmodel,'-', 'LineWidth', 3')
    plot(ytime, ydata + 1.96*sigma, 'b-')
    plot(ytime, ydata - 1.96*sigma, 'b-')
    xlabel ('time (hours)')
    ylabel(' N(t)')
    legend ('Untreated control data', 'model fit')
    legend box off
    title ('Fit for g_{s} & K using untreated control')
    set(gca,'FontSize',20,'LineWidth',1.5)
 %% Now use this and fit for gr, ds, and alpha from single pulse treatments
   
   

