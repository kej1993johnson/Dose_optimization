% This script tests the calibration method by randomly sampling parameter
% space, generating in silico data, and running the calibration for 1. N(t)
% only, 2. N(t) and phi(t) and 3. N(t) and only 3 phi(t) time points.


% It uses the N0s, time vectors, and error in data for weighting from the
% actual data.
 close all; clear all; clc
 
 %% Load in data
S = load('../out/trajfit.mat');
traj= S.traj;
% Separate by dose
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
        date = {'8-16-18'}; % pull from the same experiment: first treat
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
    Nfin = 5.5e4;
    else
    Nfin = 4e4;
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
kdrug = 0.0175;
k = 0.5;
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
%% Set & store known parameters
nsamps = 100;

pbounds = [0, 0.01; 0, 1; 0,1]; 
phidom =linspace(1,1,100);
rsdom = linspace(0.001, 0.1, 100);
carcapdom = linspace(4.8e4, 5.5e4, 100);
alphadom = linspace(0, 0.01, 100);
rrdom = linspace(0, 0.001, 100);
dsdom = linspace(0,0.01, 100);
%%
for i = 1:nsamps
    phi0=randsample(phidom,1);
    rs = randsample(rsdom,1);
    carcap = randsample(carcapdom,1);
    alpha = randsample(alphadom,1);
    rr = randsample(rrdom,1);
    ds = randsample(dsdom,1);
    dr = 0;
    pallstore(i,:) = [phi0, rs,carcap, alpha, rr, ds ,dr];
    puntfstore(i,:) = [ rs, carcap];
    pfitstore(i,:) = [ alpha, rr, ds];
end
%% Generate untreated control data
sigmaunt = trajsum(1).Nstd(1:end);
ytimeunt = trajsum(1).tvec;
Uunt = trajsum(1).U;
N0unt = trajsum(1).Nmean(1);
tdrug = 1;
eta=500;

for i = 1:nsamps
    P= num2cell(pallstore(i,:));
    [phi0, rs, carcap, alpha, rr, ds ,dr]= deal(P{:});
    S0 = N0unt*phi0;
    R0 = N0unt*(1-phi0);
    P= num2cell(pallstore(i,:));
    pit = [S0,R0, rs, carcap, alpha, rr, ds, dr];
    % Generate in silico data and store it
    [Nsrunt, ~,~] = fwd_Greene_model(pit, ytimeunt, Uunt, dt, tdrug);
    Nunt = Nsrunt(:,1) + normrnd(0, eta,[length(ytimeunt) 1]);
    phiunt = Nsrunt(:,2)./Nunt;
    Nstoreunt(:,i)= Nunt;
    phistoreunt(:,i) = phiunt;
    
    % Fit to insilico data and store fit parameters
    [punt, Nmodunt] = fit_untreated(Nunt,ytimeunt, sigmaunt);
    puntfit(i,:) = punt;
    Nfitunt(:,i) = Nmodunt;
    CCC_N(i)= f_CCC([Nunt, Nmodunt], 0.05);
    CCC_punt(i) = f_CCC([puntfstore(i,:)', punt'], 0.05);
end   
figure;
ind = 1;
plot(ytimeunt, Nstoreunt(:,ind), '*')
hold on
plot(ytimeunt, Nfitunt(:,ind), '-')
plot(ytimeunt, Nstoreunt(:,ind) + 1.96*sigmaunt, 'k-')
plot(ytimeunt, Nstoreunt(:,ind) - 1.96*sigmaunt, 'k-')
text(ytimeunt(40), Nfitunt(30,ind), ['CCC_{N}=', num2str(CCC_N(ind)),', CCC_{punt}=', num2str(CCC_punt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
legend ('in silico data', 'model fit', '95% CI on data', 'Location', 'NorthWest')
legend box off
title ('Example fit to untreated control')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Generate dosed data and fit it using puntfit

psetID = [1, 2, 3, 7];
pfitID = [4, 5, 6];
% Get what we need from real data
sigmafit = [];
ytimefit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
for i = 2:6%length(trajsum)
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);

% Generate the dosed data 
for i = 1:nsamps
    
    P=num2cell(pallstore(i,:));
    [phi0, rs, carcap, alpha, rr, ds, dr]= deal(P{:});
    pset = [phi0, rs, carcap, dr];
    pfset = [alpha, rr, ds];
    Ntrt = [];
    phitrt = [];
    for j = 2:6
        Nsr = [];
        U = trajsum(j).U;
        tvec = trajsum(j).tvec;
        N0 = trajsum(j).Nmean(1);
        S0 = phi0*N0;
        R0 = (1-phi0)*N0;
        pit = [S0, R0, rs, carcap, alpha, rr, ds, dr]; 
        [Nsr, ~, ~] = fwd_Greene_model(pit, tvec, U, dt, tdrug);
        Nsr = Nsr + normrnd(0, eta,[length(tvec) 3]);
        Ntrt = vertcat(Ntrt, Nsr(:,1));
        phitrt = vertcat(phitrt, Nsr(:,2)./Nsr(:,1));
    end
    Ntrtstore(:,i) = Ntrt;
    phitrtstore(:,i) = phitrt;
    
    % Now fit your in silico data
    % first fit Ntrt
    rsfit = puntfit(i,1);
    carcapfit = puntfit(i,2);
    pset = [phi0, rsfit, carcapfit, dr];
    theta = [0.0035, 0.1*rsfit, 0.001];
    [pbestf,N_model, negLL] = fit_fxn_Greene(Ntrt,sigmafit, pfitID, psetID, theta, pset, ytimefit, Uvec, lengthvec, N0s, pbounds);
    pfittrt(i,:) = pbestf;
    Nfittrt(:,i) = N_model;
    CCC_Ntrt(i) = f_CCC([N_model, Ntrt], 0.05);
    CCC_ptrt(i) = f_CCC([pfitstore(i,:)', pbestf'], 0.05);
    

end
%%
ind = 80
figure;
plot(ytimefit, Ntrtstore(:,ind), '*', 'LineWidth',2)
hold on
plot(ytimefit, Nfittrt(:,ind), '.', 'LineWidth',2)
plot(ytimefit, Ntrtstore(:,ind)+1.96*sigmafit, 'k.')
plot(ytimefit, Ntrtstore(:,ind)-1.96*sigmafit, 'k.')
text(ytimeunt(20), Nfitunt(20,ind), ['CCC_{Ntrt}=', num2str(CCC_Ntrt(ind)),', CCC_{pfit}=', num2str(CCC_ptrt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
legend ('in silico data', 'model fit', '95% CI on data', 'Location', 'NorthWest')
legend box off
title ('Example of fit to simulated treated data')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Find accuracy metrics

CCC_Ntrtmean = mean(CCC_Ntrt)
CCC_Nuntmean=mean(CCC_N)
CCC_puntmean = mean(CCC_punt)
CCC_ptrtmean = mean(CCC_ptrt)

% Find average percent parameter error
for i = 1:nsamps
    param_err_trt(i,:) = abs(pfitstore(i,:)-pfittrt(i,:))./pfitstore(i,:);
    param_err_unt(i,:) = abs(puntfstore(i,:) - puntfit(i,:))./puntfstore(i,:);
end
for j = 1:3
index = param_err_trt(:,j) == Inf
avg_param_err_trt(j) = median(param_err_trt(index ==0,j));
end
std_param_err_trt=std(param_err_trt)
avg_param_err_unt = nanmean(param_err_unt)
std_param_err_unt=std(param_err_unt)


avg_param_err_all = mean(avg_param_err_trt)
avg_punt_err_all = mean(avg_param_err_unt)
figure;
subplot(1, 3, 1)
plot(1:1:nsamps, pfitstore(:,1))
hold on
plot(1:1:nsamps, pfittrt(:,1))
xlabel ('Simulation Run')
ylabel('\alpha value')
legend('true', 'estimated')
title(['\alpha ', num2str(round(avg_param_err_trt(1)*100,1)),' % error'])
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1, 3, 2)
plot(1:1:nsamps, pfitstore(:,2))
hold on
plot(1:1:nsamps, pfittrt(:,2))
xlabel ('Simulation Run')
ylabel('rr value')
legend('true', 'estimated')
title(['rr ',num2str(round(avg_param_err_trt(2)*100,1)),' % error'])
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1, 3, 3)
plot(1:1:nsamps, pfitstore(:,3))
hold on
plot(1:1:nsamps, pfittrt(:,3))
xlabel ('Simulation Run')
ylabel('ds value')
legend('true', 'estimated')
title(['ds ',num2str(round(avg_param_err_trt(3)*100,1)),' % error'])
set(gca,'FontSize',20,'LineWidth',1.5)