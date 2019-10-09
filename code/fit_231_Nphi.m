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
[rsf, carcapNf, alphaf, rrf, dsf, drf, k, kdrugi, gtot] = deal(P{:});
% We will use this carcapNf only, and the k and kdrug to be consistent
kdrug = kdrugi*0.75; % make the drug last longer
% load in the scRNAseq estimates of phi(t) and the corresponding time
% vectors from python

phi_est_filename = '../data/phi_t_est_pyth.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t;
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells;
sigtech = 1e-2;
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

%% Plot the phi(t) data from single cell sequencing output
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
plot(tgen, Ub, 'b-', 'LineWidth', 3)
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
 % Again, clean data for fitting... (cut off data before carcap)
 for i = 1:length(trajsum)
%     if i ==1
%     Nfin = 5e4;
%     else
    Nfin = 3.5e4;
    %end
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
Cdoxvec = [];
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    Cdoxvec = vertcat(Cdoxvec,Cdox);
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
  dosevec = [ 1, 4, 7, 9];
for j=1:length(dosevec)
    i = dosevec(j);
 
     subplot(2,1,2)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 16)
         plot(trajsum(i).tvec, trajsum(i).Nmean + 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) data')
        set(gca,'FontSize',20,'LineWidth',1.5)
        
        dt = 1;
       subplot(2,1,1)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'*', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('Effective dose (u(t)) ')
        set(gca,'FontSize',20,'LineWidth',1.5)
        legend('0 nM', '75 nM', '200 nM', '500 nM', 'Location', 'NorthEast')
        legend boxoff
save('../out/trajsumfit231.mat', 'trajsum')
end
%% Compile dosed data and fit it using puntfit
% Here we're going to generate dosed data and output N(t) and phi(t)
% Change pset to only phi0 and carcap
psetID = [ 3, 4]; % phi0, carcapN, carcapphi
pfitID = [ 1, 2, 5, 6, 7, 8]; % corresponds to phi0, rs, alpha, rr, ds, dr

% Get what we need from real data
sigmafit = [];
ytimefit = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
dosevec = [ 1, 4, 7, 9];
for j=1:length(dosevec)
    i = dosevec(j);
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);
Cdoxdata = Cdoxvec(dosevec);
%% Now fit your data using both Ntrt and phitrt

% Currently, we are going to pretend to set rr. In the future, we could
% also set or just bound ds or dr
rr_to_pop_ratio= 0.03/0.0392;
zrdata = rr_to_pop_ratio;
zddata = 0.1;
%pset = [phi0, carcapNf, carcapphi, zrdata];
pset = [ carcapNf, carcapphi];
rstar = phi0/(1-phi0);
zrguess = 0.4;
rsguess =  1.8*gtot;
alphaguess = 0.1;
dsguess = 0.04;
zdguess = 0.1;
phi0guess = phi0;
theta = [phi0guess, rsguess, zrguess, alphaguess, dsguess, zdguess];

%phi0, rs, alpha, zr, ds, and zd
% Unfortunately it appears that the 
%pbounds =  [ 0.5*gtot, 2; 0,1; 0,2; 0.0,0.1; 0,1]; 
pbounds =  [ 0,1; 0,1; 0,1; 0,1; 0,1; 0,1]; 
%Give this function both Ntrt and phitrt
% can toggle the amount that we weigh each portion..
% lambda =0-- fit on N(t) only. if lambda =1, fit on phi(t) only
lambda = 0.995;
% This function internally instead of actually fitting rr and dr, fits the ratio 
[pbest,N_model, phi_model, negLL, err_N, err_phi] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% Adjust transformed parameter estimates so we capture estimated value
% of rr and dr (since these are what we saved).
phi0 = pbest(1);
rs =pbest(2);
alpha = pbest(3);
zr = pbest(4);
rr = rs*zr;
ds = pbest(5);
zd = pbest(6);
dr = ds*zd;

% Simulate the phi(t) trajectory from all time points
phi_model_long = simmodelgreenephi2(pbest, tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);

CCC_vec(1) = f_CCC([N_model, ydatafit], 0.05);
CCC_vec(2) = f_CCC([phi_model, phitrt], 0.05);
fvalbest = negLL
pbest = pbest

psave = [phi0, carcapNf, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot];
save('../out/ptest.mat', 'psave')
%% Plot first pass fitting results with an arbitrary lambda

figure;
subplot(1,2,1)
errorbar(ytimefit, ydatafit, 1.96*sigmafit/2, 'b*')
hold on
plot(ytimefit, N_model, 'k.', 'LineWidth',1)
for j = 1:length(dosevec)
    i = dosevec(j);
 text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'])
end
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
plot(tgen, phi_model_long,'k-', 'LineWidth',2)
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

%% Model fit and data 
for i= 1:length(trajsum)
Nsri= [];
Nsrlow = [];
Nsrhigh = [];
dt = 1;
tdrug = 1;
Cdox = trajsum(i).Cdox;
tgeni = 0:1:trajsum(i).tvec(end);
Udata=k*Cdox*exp(-kdrug*(tgeni))/(0.1*Cdoxmax);
N0 = trajsum(i).Nmean(1);
p = [ phi0, carcapNf,rs,alpha, zr, ds, zd];
[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, trajsum(i).tvec, N0, Udata, dt, tdrug);

% Add to the trajsum data structure the model fit for that Cdox (whether it
% is calibrated or predicted)
trajsum(i).Nsrmod = Nsri;

end
trajsum(1).color = [0.4 0.6 0];
figure;
subplot(1,2,1)
for j = 1:length(dosevec)
    i = dosevec(j);
    plot(trajsum(i).tvec, trajsum(i).Nsrmod(:,1), 'color', trajsum(i).color, 'LineWidth', 2)
    hold on
end
for j = 1:length(dosevec)
    i = dosevec(j);
    plot(trajsum(i).tvec, trajsum(i).Nmean,'.', 'color', trajsum(i).color, 'LineWidth', 2)
    text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 16, 'color', trajsum(i).color)
    plot(trajsum(i).tvec, trajsum(i).Nmean + 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
    plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
end
xlabel ('time (hours)')
ylabel(' N(t)')
%legend ( 'N(t) model fit', 'N(t)data & 95% CI', 'Location', 'NorthWest')
%legend box off
title ('N(t) data and model calibration')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
errorbar(tbot, phitrt, phisigfit/2,  'go', 'LineWidth', 4)
plot(tgen, phi_model_long,'g-', 'LineWidth',2)
%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi_{sens}(t)')
%legend ('\phi_{sens}(t) data','model fit','95% CI on data', 'Location', 'NorthWest')
%legend box off
title ('\phi(t) data and model calibration')
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([0 1.0])
xlim([0 1656])
%% Now we want to find the pareto surface by fitting at different values of lambda
% at each lambda capture the error in N(t) and phi(t) and plot for
% different values of lambda. 

% Eventually, we want to get more of these
lambdavec1 = linspace(0.99, 1, 100);
% instead try randomly sampling lambda
lambdavec  = 0.99 + (0.01)*rand(100,1);
%%
for i = 1:length(lambdavec)
    lambdai = lambdavec(i);
    [pbesti(:,i),N_modeli(:,i), phi_modeli(:,i), negLLi(:,i), weightederr_Ni, weightederr_phii] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambdai, pbounds);
    sumsqerrN(i) =weightederr_Ni;
    sumsqerrphi(i) =weightederr_phii;
    % Also calculate CCC in N and CCC in phi
    CCC_phi(i) = f_CCC([phi_modeli(:,i), phitrt], 0.05);
    CCC_N(i) = f_CCC([N_modeli(:,i), ydatafit], 0.05);
end
%%
pareto = horzcat(sumsqerrphi', sumsqerrN', pbesti')
save('../out/pareto.mat', 'pareto')

%%
paramnames = {'\phi_{0}','rs', 'zr', '\alpha', 'ds', 'zd'};
figure;
for i = 1:length(pbest)
    subplot(1,length(pbest), i)
    plot(lambdavec, pbesti(i,:), '.')
    hold on
    xlim([min(lambdavec), max(lambdavec)])
    xlabel('\lambda')
    ylabel('parameter value')
    title(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5)
end
%% Plot error in N vs. error in phi
figure;
scatter(sumsqerrphi, sumsqerrN,[],lambdavec, 'filled')
hold on
plot(err_phi, err_N, 'r*', 'LineWidth', 4)
colorbar
xlabel('weighted sum squared error in \phi(t)')
ylabel('weighted sum squared error in N(t)')
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
xlim([ 0 1e3])
ylim([ 3e3 5e4])
legend('lambda iterations', 'current best fit')
legend boxoff
title('Objective function values')
%%
figure;
scatter(CCC_phi, CCC_N,[],lambdavec, 'filled')
hold on
plot(CCC_vec(2), CCC_vec(1), 'r*', 'LineWidth', 4)
colorbar
xlabel('CCC in \phi(t)')
ylabel('CCC in N(t)')
xlim([ 0.88 1])
ylim([0.88 1])
set(gca,'FontSize',20,'LineWidth',1.5)
legend('lambda iterations', 'current best fit')
legend boxoff
title('CCC in model vs. data')
%% Select the parameter sets that are above the CCC threshold in N and phi fit
iphi = CCC_phi>0.88;
iN = CCC_N>0.88;
ikeep = and(iphi, iN);
phi0vals = pbesti(1,ikeep)';
phi0vals = vertcat(phi0vals, phi0);
rsvals = pbesti(2,ikeep)';
rsvals = vertcat(rsvals,rs);
alphavals = pbesti(3,ikeep)';
alphavals = vertcat(alphavals,alpha);
rr_rs_ratio= pbesti(4,ikeep)';
rr_rs_ratio = vertcat(rr_rs_ratio, zr);
dsvals = pbesti(5,ikeep)';
dsvals = vertcat(dsvals, ds);
dr_ds_ratio = pbesti(6,ikeep)';
dr_ds_ratio = vertcat(dr_ds_ratio, zd);
%%
lambdas = lambdavec(ikeep);
lambdas = vertcat(lambdas, lambda);
err_Nvals = sumsqerrN(ikeep)';
err_Nvals = vertcat(err_Nvals, err_N);
err_phivals = sumsqerrphi(ikeep)';
err_phivals = vertcat(err_phivals, err_phi);
CCC_Nvals = CCC_N(ikeep)';
CCC_Nvals = vertcat(CCC_Nvals, CCC_vec(1));
CCC_phivals = CCC_phi(ikeep)';
%%
CCC_phivals = vertcat(CCC_phivals, CCC_vec(2));

% Pass these parameter sets to a table to be imported into the prediction
% script
param_table = table(phi0vals,rsvals, alphavals, rr_rs_ratio, dsvals, dr_ds_ratio, err_Nvals, err_phivals, CCC_Nvals, CCC_phivals, lambdas )
writetable(param_table,'../out/pbest_table.csv')
% test we can load it back in
pareto = readtable('../out/pbest_table.csv')
%%
paramnames = {'\phi_{0}','rs', 'zr', '\alpha', 'ds', 'zd'};
figure;
for i = 1:length(pbest)
    subplot(2,length(pbest), i)
    scatter(param_table.lambdas, param_table{:, i},[],param_table.CCC_Nvals, 'filled')
    colorbar
    hold on
    xlim([min(lambdas), max(lambdas)])
    xlabel('\lambda')
    ylabel('parameter value')
    title(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(2,length(pbest), i+length(pbest))
    scatter(param_table.lambdas, param_table{:, i},[],param_table.CCC_phivals, 'filled')
    colorbar
    hold on
    xlim([min(lambdas), max(lambdas)])
    xlabel('\lambda')
    ylabel('parameter value')
    title(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5)
end

%%
% Next up-- how do we find the optimal value along this curve? And once we
% do that, how do we set lambda accordingly?
% Probably need to do some literature digging for this...

figure;
for i = 1:length(pbest)
    subplot(2,length(pbest), i)
    plot(pbesti(i,:), sumsqerrN, 'g.')
    set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
    %xlim([min(paretomat(:,i)) max(paretomat(:,i))])
    xlabel([paramnames(i)])
    ylabel('error in N')
    
    subplot(2,length(pbest), i+length(pbest))
    plot(pbesti(i,:), sumsqerrphi, 'r.')
    set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
    %xlim([min(paretomat(:,i)) max(paretomat(:,i))])
    xlabel([paramnames(i)])
    ylabel('error in \phi')
    
end

%% Make a figure of parameter relationship
np = length(pbest);
itvec = 1:1:np*np;
ct_iter = 0;
figure;
for i = 1: length(pbest)
    for j = 1:length(pbest)
        ct_iter = ct_iter+1;
        subplot(length(pbest), length(pbest),itvec(ct_iter))
        scatter(pbesti(i,:), pbesti(j,:),10, sumsqerrN, 'filled')
        colorbar
        set(gca,'FontSize',10,'LineWidth',1.5)
        xlabel([pbestnames(i)])
        ylabel([pbestnames(j)])
        if i ==j
            hist(pbesti(i,:))
            xlabel([pbestnames(i)])
        end
    end
end
%%
ct_iter = 0;
figure;
for i = 1: length(pbest)
    for j = 1:length(pbest)
        ct_iter = ct_iter+1;
        subplot(length(pbest), length(pbest),itvec(ct_iter))
        scatter(pbesti(i,:), pbesti(j,:),10, sumsqerrphi, 'filled')
        colorbar
        set(gca,'FontSize',10,'LineWidth',1.5)
        xlabel([pbestnames(i)])
        ylabel([pbestnames(j)])
        if i ==j
            hist(pbesti(i,:))
            xlabel([pbestnames(i)])
        end
    end
end


%%
ct_iter = 0;
figure;
for i = 1: length(pbest)
    for j = 1:length(pbest)
        ct_iter = ct_iter+1;
        subplot(length(pbest), length(pbest),itvec(ct_iter))
        scatter(pbesti(i,:), pbesti(j,:),10, negLLi, 'filled')
        colorbar
        set(gca,'FontSize',10,'LineWidth',1.5)
        xlabel([pbestnames(i)])
        ylabel([pbestnames(j)])
        if i ==j
            hist(pbesti(i,:))
            xlabel([pbestnames(i)])
        end
    end
end

%% Profile likelihood analysis


% We are going to write a function that performs the profile likelihood
% analysis by iterating through the pbest, setting that value (using your
% psetID and changing pset), letting all the other parameters float, and
% finding the best fitting parameters and the 