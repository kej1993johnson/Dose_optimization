% First pass at fitting the N(t) data and the phi(t) estimates (just
% transcribed in from the latest python results from scRNAseq)



% It uses the N0s, time vectors, and error in data for weighting from the
% actual data.
 close all; clear all; clc
 
 %% Load in data from 231 dose response and from fit
S = load('../out/trajfit231.mat');
traj= S.traj;

% Use the fit parameters from N(t) only to make the parameter domains
p4fit = load('../out/p4fit231.mat');
p4fit = struct2cell(p4fit);
p4fit = cell2mat(p4fit);
P = num2cell(p4fit);
% Can use some of these as first guesses/ballpark ranges of what to expect
[rsf, carcapNf, alphaf, rrf, dsf, drf, k, kdrug] = deal(P{:});

% load in the scRNAseq estimates of phi(t) and the corresponding time
% vectors from python
phi_est_filename = '../data/phi_t_est.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t;
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells
sigtech = 1e-4;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
N0phi = 0.8*0.24e6; % set this because this is what we think we seeded for this experiment
phi0= phitrt(1);
S0phi = phi0*N0phi;
R0phi = (1-phi0)*N0phi;
Cdoxphi = 550;
Cdoxmax = 1000;
tgen = [0:1:tbot(end)];
Ub=k*Cdoxphi*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
lengthvecphi = [length(tbot), length(tgen)];
carcapphi =20e6; % set this based on final outgrowth 

%% Plot the phi(t) data from single cell
figure;
subplot(2,1,1)
errorbar(tbot, phitrt, phisigfit/2,  'go', 'LineWidth', 4)
hold on
errorbar(tbot, 1-phitrt, phisigfit/2, 'ro', 'LineWidth', 4)
legend('\phi_{sens}(t)', '\phi_{res}(t)', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel(' \phi(t) data')
xlim([ 0 1656])
ylim([ 0 1.2])
title('\phi(t) for dosing for scRNAseq expt')
subplot(2,1,2)
plot(tgen, Ub, 'k-', 'LineWidth', 3)
xlabel('time (hours)')
ylabel('Effective dose U(t)')
title('U(t) for dose for scRNAseq expt')
xlim([ 0 1656])


%% Make cleaned N(t) data from range of dox concentrations
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
% Make a new structure which combines each dox concentration
 for i = 1:length(uniqdose)
    trajsum(i).Cdox = [];
    trajsum(i).Nmat = [];
    trajsum(i).nreps = 0;
    trajsum(i).tmat = [];
 end

 for i = 1:length(uniqdose) % number of unique seed numbers
    for j = 1:length(traj)
        date = {'5-6-19'}; % pull from the same experiment: first treat
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
 % Again, clean data for fitting...
 for i = 1:length(trajsum)
    if i ==1
    Nfin = 5e4;
    else
    Nfin = 3.5e4;
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
% Test and set U(t) curves
dt = 1;
% input time vectors for each different dose response
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest))/(0.1*Cdoxmax);  
    subplot(1, length(trajsum),i)
    plot(ttest, trajsum(i).U, 'b-', 'LineWidth',2)
    ylim([0 5])
    xlim([ 0 ttest(end)])
    xlabel('time (hours)')
    ylabel('Effective dose (U(t))')
    title([num2str(trajsum(i).Cdox),' nM Dox'])
end
% Now find mean and standard deviation vectors
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
         plot(trajsum(i).tvec, trajsum(i).Nmean + 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
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

%% Fit the untreated control data
% Fit this only to N(t) 
sigmaunt = trajsum(1).Nstd(1:end);
ytimeunt = trajsum(1).tvec;
Uunt = trajsum(1).U;
N0unt = trajsum(1).Nmean(1);
tdrug = 1;
Nunt = trajsum(1).Nmean;

S0 = N0unt*phi0;
R0 = N0unt*(1-phi0);
% Fit to insilico data and store fit parameters
[punt, Nmodunt] = fit_untreated(Nunt,ytimeunt, sigmaunt);
% punt is gtot, carcap1
% gtot is not a parameter we set, so dont include this when comparing
% parameter error but fine to keep track of it.
gtot = punt(1);
carcapN = punt(2);

% Plot the fit
figure;
plot(ytimeunt, Nunt, '*')
hold on
plot(ytimeunt, Nmodunt, '-')
plot(ytimeunt, Nunt + 1.96*sigmaunt, 'k-')
plot(ytimeunt, Nunt - 1.96*sigmaunt, 'k-')
%text(ytimeunt(20), Nfitunt(15,ind), ['CCC_{N}=', num2str(round(CCC_N(ind),3)),', % error_{carcap}=', num2str(round(pct_err_carcap(ind),2)), '%'],'FontSize',14)
xlabel ('time (hours)')
ylabel(' N(t)')
legend ('untreated data', 'model fit', '95% CI on data', 'Location', 'NorthWest')
legend box off
title ('Fit to untreated control  data')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Compile dosed data and fit it using puntfit
% Here we're going to generate dosed data and output N(t) and phi(t)
% Change pset to only phi0 and carcap
psetID = [1, 3, 4]; % phi0, carcapN, carcapphi
pfitID = [2, 5, 6, 7, 8]; % corresponds to rs, alpha, rr, ds, dr
% Get what we need from real data
sigmafit = [];
ytimefit = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
for i = 2:6%length(trajsum)
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);
%% Now fit your in  data using both Ntrt and phitrt
pset = [phi0, carcapN, carcapphi];
rstar = phi0/(1-phi0);
zrguess = 0.2;
rsguess =  gtot;
alphaguess = 0.015;
dsguess = 0.015;
zdguess = 0.1;
theta = [rsguess, alphaguess, zrguess, dsguess, zdguess];
%rs, zr, alpha, ds, and zd
pbounds = [0,5*gtot; 0,1; 0,1; 0,1; 0,1]; 
%Give this function both Ntrt and phitrt
% can toggle the amount that we weigh each portion..
% lambda =0-- fit on N(t) only. if lambda =1, fit on phi(t) only
lambda = 0.9;
% This function internally instead of actually fitting rr and dr, fits the ratio 
[pbest,N_model, phi_model, negLL, err_N, err_phi] = fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% Adjust transformed parameter estimates so we capture estimated value
% of rr and dr (since these are what we saved). 
rr = pbest(1)*pbest(3) %rr= zr*rs
dr = pbest(4)*pbest(5) %dr = zd*dr
rs = pbest(1)
alpha = pbest(2)
ds = pbest(4)
% Simulate the phi(t) trajectory from all time points
phi_model_long = simmodelgreenephi2(pbest, tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);

CCC_vec(1) = f_CCC([N_model, ydatafit], 0.05)
CCC_vec(2) = f_CCC([phi_model, phitrt], 0.05)

%% Plot first pass fitting results with an arbitrary lambda

figure;
subplot(1,2,1)
errorbar(ytimefit, ydatafit, 1.96*sigmafit/2, 'b*')
hold on
plot(ytimefit, N_model, 'k.', 'LineWidth',1)
%plot(ytimefit, ydatafit-1.96*sigmafit, 'b.')
%text(ytimeunt(20), Nfitunt(20,ind), ['CCC_{Ntrt}=', num2str(CCC_Ntrt(ind)),', CCC_{pfit}=', num2str(CCC_ptrt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
legend ( 'N(t) data & 95% CI','model fit', 'Location', 'NorthWest')
legend box off
title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
errorbar(tbot, phitrt, phisigfit/2,  'go', 'LineWidth', 4)
plot(tgen, phi_model_long,'k-', 'LineWidth',1)
%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi_{sens}(t)')
legend ('\phi_{sens}(t) data','model fit','95% CI on data', 'Location', 'NorthWest')
legend box off
title (['\phi(t), CCC_{\phi}=', num2str(CCC_vec(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([0 1.5])
xlim([0 1656])
%pause
%% Now we want to find the pareto surface by fitting at different values of lambda
% at each lambda capture the error in N(t) and phi(t) and plot for
% different values of lambda. 

lambdavec = linspace(1e-4, 1-(1e-4), 100);

for i = 1:length(lambdavec)
    lambdai = lambdavec(i);
    [pbesti(:,i),N_modeli(:,i), phi_modeli(:,i), negLLi(:,i), err_Ni, err_phii] = fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambdai, pbounds);
    sumsqerrN(i) =sum(err_Ni.^2);
    sumsqerrphi(i) = sum(err_phii.^2);
end
%% Plot error in N vs. error in phi
figure;
plot(sumsqerrphi, sumsqerrN, 'ro', 'LineWidth', 2)
xlabel('sum squared error in \phi(t)')
ylabel('sum squared error in N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)%, 'Xscale', 'log', 'Yscale', 'log')
xlim([ 1e-9 2e-4])
ylim([ 4e9 9e9])
title('Pareto Optimality')

% Next up-- how do we find the optimal value along this curve? And once we
% do that, how do we set lambda accordingly?
% Probably need to do some literature digging for this...
%% Profile likelihood analysis


% We are going to write a function that performs the profile likelihood
% analysis by iterating through the pbest, setting that value (using your
% psetID and changing pset), letting all the other parameters float, and
% finding the best fitting parameters and the 