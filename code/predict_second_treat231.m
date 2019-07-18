% This script loads in the 231 single pulse treatment parameter estimates
% and makes a prediction for a second dose treatment that we actually have
% data on...

% We start by simulating completely the effect of a second treatment using
% the forward model and the fitted parameters. We choose a dosing schedule
% that has been observed in our data.

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/trajfit231.mat');
traj= S.traj;

S = load('../out/trajsumfit231.mat');
trajsum = S.trajsum;

p4fit = load('../out/p4fit231.mat');
p4fit = struct2cell(p4fit);
p4fit = cell2mat(p4fit);
P = num2cell(p4fit);
[rs, carcap, alpha, rr, ds, dr] = deal(P{:});
%% Run this forward for a single pulse treatment at 75 nM 
% Again, assume R0=0; dr = 0 and
kdrug = 0.0175;
k = 0.7;
tdrug =1; % stands for first treatment
dt = 1;
phi0=0.57; % from scRNASeq pre-treatment

% Pull from trajsum dose = 75
for i = 1:length(trajsum)
    if trajsum(i).Cdox == 75
        tvec = trajsum(i).tvec;
        Nmean = trajsum(i).Nmean;
        Nstd = trajsum(i).Nstd;
        U1 = trajsum(i).U;
        S0 = phi0*trajsum(i).Nmean(1);
    end

end
tvecdat = tvec;
R0 = (1-phi0)*trajsum(i).Nmean(1);
p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
tlong = 0:1:1344;
tvec = tlong; % simulate long term treatment dynamics
Cdox = 75;
U1 = k*Cdox*exp(-kdrug*(tlong));
[Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U1, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);
figure;
plot(tvecdat, Nmean,'k*', 'LineWidth', 2)
hold on
plot(tvec, Nsri(:,1), 'LineWidth',2)
plot(tvecdat, Nmean + Nstd, 'color', 'k')
plot(tvecdat, Nmean - Nstd, 'color', 'k')
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) for single 75 nM pulse treatment')
legend ('data mean', 'model', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tvecdat, Nmean,'k*', 'LineWidth', 2)
hold on
plot(tvec, Nsri(:,1), 'LineWidth',2)
plot(tvecdat, Nmean + Nstd, 'color', 'k')
plot(tvecdat, Nmean - Nstd, 'color', 'k')
plot(tvec, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tvec, Nsri(:,2), 'g','LineWidth',2)
plot(tvec, Nsri(:,3),'r', 'LineWidth',2)
text(tvec(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)])
text(tvec(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 75 nM pulse treatment')
legend ('N data', 'Upper bound', 'Lower bound','N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tvec, sensfrac_t,'g', 'LineWidth',3)
hold on
plot(tvec, resfrac_t, 'r','LineWidth',3)
text(tvec(end-100), sensfrac_t(end-100), ['\phi_{sens_t=8WPT}=', num2str(fracSens)], 'FontSize', 14)
text(tvec(end-100), resfrac_t(end-100), ['\phi_{res_t=8WPT}=', num2str(1-fracSens)], 'FontSize', 14)
plot(tvec(end-100), sensfrac_t(end-100), 'k*', 'LineWidth',5)
plot(tvec(end-100), resfrac_t(end-100), 'k*', 'LineWidth',5)
text(tvec(30), sensfrac_t(30), ['\phi_{sens_t=30hrs}=', num2str(sensfrac_t(30))], 'FontSize', 14)
text(tvec(30), resfrac_t(30), ['\phi_{res_t=30hrs}=', num2str(resfrac_t(30))], 'FontSize', 14)
plot(tvec(30), sensfrac_t(30), 'k*', 'LineWidth',5)
plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
text(tvec(504), sensfrac_t(504), ['\phi_{sens_t=3WPT}=', num2str(sensfrac_t(504))], 'FontSize', 14)
text(tvec(504), resfrac_t(504), ['\phi_{res_t=3WPT}=', num2str(resfrac_t(504))], 'FontSize', 14)
plot(tvec(504), sensfrac_t(504), 'k*', 'LineWidth',5)
plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
xlim([ 0 tlong(end)])
ylabel ('Proportion of cells')
title('\phi_{S} & \phi_{R} following single 75 nM pulse treatment')
% legend ( '\phi_{S}','\phi_{R}', 'Location', 'NorthWest')
% legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

%% Test function to create an easy to work with combined data structure
% Note for now this isn't going to work since we don't have the 231 repeat
% treatment data
filter_criteria = 'treatmentnum';
ntreat = 2;
dose = 75;
date = '12-19-18';
[traj2] = comb_data_fxn(traj, filter_criteria, ntreat, dose, date)

%% Write a function to find phi at a certain time
% If repeat dose is within the time frame measured in the first pusle
% treatment
repeat_trts = 1:1:8;
  Cdox = 150;
figure;
for i = 1:length(repeat_trts)
   % Find the Number of sensitive and resisitant cells predicted at start
   % of second treatment
    tdose = repeat_trts(i)*24*7;
    tfirst = 0:1:tdose;
    U1 = k*Cdox*exp(-kdrug*(tfirst));
    p1 = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p1, tfirst, U1, dt, tdrug);
    ind = find(ismember(tfirst, tdose), 1, 'first');
    fracS(i)= Nsri(ind,2)/Nsri(ind,1);
    % Set the initial number of the total cells and the sensitive and
    % resistant cells at the time of repeat treatment
    N01(1,i) = 2e3; % lets say that we simulated reseeding 2000 cells. 
    % We assume the proportion of sensitive and resistant cells remains the
    % same
    S01(1,i) = fracS(i)*N01(1,i);
    R01(1,i) = (1-fracS(i))*N01(1,i);


% Simulate second dose 
    tfin = 672; % this is how long we want to monitor the repeat dose (about 5 weeks)
    tvec2 = 0:dt:tfin;
    U2 = k*Cdox*exp(-kdrug*(tvec2));
    tdrug = 1; % only simulating one dose)
    p2 = [S01(1,i), R01(1,i), p(3:end)];
    [Nsr2, tcrit2, Ncrit2] = fwd_Greene_model(p2, tvec2, U2, dt, tdrug);
    Nmodpred = Nsr2;
    tvecmod = tvec2;




subplot(2, length(repeat_trts)/2, i)

    hold on
    plot (tvecmod, Nmodpred(:,1),'b-' ,'LineWidth',2)
    plot(tvecmod, Nmodpred(:,2), 'g-', 'LineWidth',2)
    plot(tvecmod, Nmodpred(:,3), 'r-', 'LineWidth', 2)
    title (['WPT=', num2str(repeat_trts(i))])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    title(['WPT=', num2str(repeat_trts(i))])
    legend ('N(t)', 'S(t)', 'R(t)', 'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    xlim([0 tfin])
    ylim([0 4e4])
end
        
