% This script loads in the traj structure which contains the dose response
% for all conditions. It creates the trajsum structure which will contain a
% mean vector and U(t) vector for each different concentration tested. This
% will be used to calibrate to the model.
 close all; clear all; clc
 
 %% Load in data structure 
S = load('../out/trajfit.mat');
traj= S.traj;

%% Separate by dose
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

 for i = 1:length(trajsum)
    if i ==1
    Nfin = 5.5e4;
    else
    Nfin = 4.5e4;
    end
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
    
    trajsum(i).tfit =round(tfit,0);
    trajsum(i).Nfit =Nfit;
    trajsum(i).tvec = trajsum(i).tfit(:,1);
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
 %% Test and set U(t) curves
% Set some arbitrary things as inputs

kdrug = 0.0175;
k = 1;
dt = 1;
% input time vectors for each different dose response
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest)); 
    subplot(1, length(trajsum),i)
    plot(ttest, trajsum(i).U, 'b-', 'LineWidth',2)
    ylim([0 300])
    xlim([ 0 ttest(end)])
    xlabel('time (hours)')
    ylabel('Effective dose (U(t))')
    title([num2str(trajsum(i).Cdox),' nM Dox'])
end
%% Now find mean and standard deviation vectors

for i = 1:length(trajsum)
    trajsum(i).Nmean = mean(trajsum(i).Nfit,2);
    trajsum(i).tvec = round(trajsum(i).tfit(:,1),0);
    trajsum(i).Nstd = std(trajsum(i).Nfit,0,2);
end
% Plot the average data 
  figure;
 for i = 1:6%length(trajsum)
     subplot(2,1,1)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'])
         plot(trajsum(i).tvec, trajsum(i).Nmean + trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - trajsum(i).Nstd, 'color', trajsum(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for different single pulse treatments')
        dt = 1;
       subplot(2,1,2)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('U(t) for different single pulse treatments')
 end


%%
% set parameters of forward model
carcap = 5.2409e4;
S0=trajsum(1).Nmean(1); % set initial conditions (assume all cells sensitive to start)
R0 = 0; 
rs = 0.0287;
rr = 0.001;
ds = 0.002;
dr = 0;
alpha =0; %0.0001;


p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];


Nsr = [];
tsum = [];
ttestsum = [];
tdrug = 1; % since we are monitoring first treatment only
ttest = [];
tvec = [];

for i = 1:6
    tvec = round(trajsum(i).tvec,0);
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    Nsr = vertcat(Nsr, Nsri);
end



figure;
plot(tsum, Nsr(:,1), 'b*')
hold on
plot(trajsum(1).tvec, trajsum(1).Nmean, 'r*')
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
N0 = ydataf(1);
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
    legend ('Untreated control data', 'model fit', 'Location', 'NorthWest')
    legend box off
    title ('Fit for g_{s} & K using untreated control')
    set(gca,'FontSize',20,'LineWidth',1.5)
 %% Now use this and fit for gr, ds, and alpha from single pulse treatments
rs = gs;
dr = 0;
carcap = carcap;
props = 1;
pset = [rs, carcap, props, dr];


%% Bayesian fit for alpha, rr, and ds using all treatments

% The simplified  model looks like this
%
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R ;
% 
% We will fit N= S +R trajectories for doses 10, 20, 35, 50, & 75
% Fit for alpha, rr, & ds


% Define transforms 
% for 3 variables
pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


% Make vectors of inputs from 2:6 of traj sum corresponding to doses desired to
% fit
sigmafit = [];
ydatafit = [];
ytimefit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];

% fit on doses 10, 20, 35, & 50 nM dox
for i = 2:5
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(2:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(2:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(2:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean)-1);
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);


% Set up forward models, fit all three nested versions of model

modelfun = @(p)simmodelgreene(p, ytimefit, N0s, pset, Uvec, lengthvec); % single exponential model with death  
%% test
pftest = [0, 0, 0*ds];
test = modelfun(pftest)
figure;
plot(ytimefit, test, '*')
hold on
plot(ytimefit, ydatafit, '*')

figure;
plot(1:1:length(Uvec), Uvec, '*')


%% INITIAL GUESSES FOR alpha, rr & ds

alphaguess =  1e-4;
rrguess = 0.015;
dsguess = 1e-3;

    % Initial guess vector
    theta = [alphaguess, rrguess, dsguess]; 

    
    % Write log likelihood function based on assumption of normally
    % distributed sampling error
    
    % Goal: maximize the probability of the data given the model. Normpdf
    % will output a probability of the data (x- 1st argument), given the
    % mean(expectation, here the model), and the variance, at each time point. 
    % take the log and minimize the NLL to maximize likeihood of data given
    % the model
    
    % Greene model
    loglikelihood = @(phat)sum(log(normpdf(yfxform(ydatafit),yfxform(modelfun(pbxform(phat))), sigmafit)));
  
   
    % Write objective functions for each model
    objfun = @(phat)-loglikelihood(phat); 
    phatbest = fminsearch(objfun, pfxform(theta));
    
    pbest = pbxform(phatbest);

    Greenemodel = modelfun(pbest);
    figure;
    plot(ytimefit, ydatafit, 'b*', 'LineWidth', 3)
    hold on
    plot(ytimefit, Greenemodel,'r*', 'LineWidth', 3')
    plot(ytimefit, ydatafit + 1.96*sigmafit, 'g.')
    plot(ytimefit, ydatafit - 1.96*sigmafit, 'g.')
    xlabel ('time (hours)')
    ylabel(' N(t)')
    legend ('data from mult treatments', 'model fit')
    legend box off
    title ('Fit for \alpha, rr and ds')
    set(gca,'FontSize',20,'LineWidth',1.5)
%%
P = num2cell(pbest); 
[alpha, rr, ds] = deal(P{:}); % our parameters

p= [N0, 0, rs, carcap, alpha, rr, ds, dr];
 figure;
 for i = 2:5%length(trajsum)
     subplot(2,1,1)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         plot(trajsum(i).tvec, trajsum(i).Nmean + trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - trajsum(i).Nstd, 'color', trajsum(i).color)
       
         xlabel('time (hours)')
        ylabel('N(t)')
        title('Model calibration to pulsed treatments 10-50 nM')
        dt = 1;
        tvec = [];
        Nsri = [];
        tvec = trajsum(i).tvec;
        U = trajsum(i).U;
         pi = p;
         pi(1) = trajsum(i).Nmean(1);
        [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
        plot(tvec, Nsri(:,1), 'color','r', 'LineWidth',2)
        plot(tcrit, Ncrit, '*', 'LineWidth',2)
        text(tcrit +5, Ncrit + 5, ['t_{crit}=', num2str(tcrit),' hours'])
        text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
        set(gca,'FontSize',20,'LineWidth',1.5)
        trajsum(i).Nmod1pulse = Nsri;



        subplot(2,1,2)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('Effective dose of each pulse treatment')
        set(gca,'FontSize',20,'LineWidth',1.5)
 end
%% Save trajsum

%p= [N0, 0, rs, carcap, alpha, rr, ds, dr];
for i = 1:length(trajsum)
    pi = p;
    pi(1) = trajsum(i).Nmean(1);
    tvec = [];
    Nsri = [];
    tvec = trajsum(i).tvec;
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
   trajsum(i).Nmod1pulse = Nsri;
    trajsum(i).params = pi;
    trajsum(i).N0 = trajsum(i).Nmean(1);
    trajsum(i).rs= rs;
    trajsum(i).carcap = carcap;
    trajsum(i).alpha = alpha;
    trajsum(i).rr = rr;
    trajsum(i).ds = ds;
    trajsum(i).dr = 0;
    trajsum(i).kdrug =kdrug;
end

    %%
p = [pset, pbest];
save('../out/pall.mat', 'p')
save('../out/trajsumfit.mat', 'trajsum')