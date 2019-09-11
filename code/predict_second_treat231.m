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

ptest = load('../out/ptest.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);

[phi0, carcapNf, carcapphi, rs, alpha, zrdata, ds, zd, k, kdrug, gtot] = deal(P{:});

%% Run this forward for a single pulse treatment at 75 nM 
% Again, assume R0=0; dr = 0 and


% Pull from trajsum dose = 75
for i = 1:length(trajsum)
    if trajsum(i).Cdox == 200
        tvec = trajsum(i).tvec;
        Nmean = trajsum(i).Nmean;
        Nstd = trajsum(i).Nstd;
        N0 = Nmean(1);
    end
end

p = [ phi0, carcapNf,rs,alpha, zrdata, ds, zd];
tlong = 0:4:1344;
dt = 1;
tdrug = 1;
Cdoxmax = 1000;
Cdox = 200;
tgen = 0:1:tlong(end);
U1=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
%tvec = tlong; % simulate long term treatment dynamics
[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, tlong, N0, U1, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);
figure;
plot(tvec, Nmean,'k*', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,1), 'LineWidth',2)
plot(tvec, Nmean + Nstd, 'color', 'k')
plot(tvec, Nmean - Nstd, 'color', 'k')
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) for single 150 nM pulse treatment')
legend ('data mean', 'model', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tvec, Nmean,'k*', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,1), 'LineWidth',2)
plot(tvec, Nmean + Nstd, 'color', 'k')
plot(tvec, Nmean - Nstd, 'color', 'k')
plot(tlong, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tlong, Nsri(:,2), 'g','LineWidth',2)
plot(tlong, Nsri(:,3),'r', 'LineWidth',2)
text(tlong(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)])
text(tlong(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 150 nM pulse treatment')
legend ('N data', 'Upper bound', 'Lower bound','N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tlong, sensfrac_t,'g', 'LineWidth',3)
hold on
plot(tlong, resfrac_t, 'r','LineWidth',3)
text(tlong(end-100), sensfrac_t(end-100), ['\phi_{sens_t=8WPT}=', num2str(fracSens)], 'FontSize', 14)
text(tlong(end-100), resfrac_t(end-100), ['\phi_{res_t=8WPT}=', num2str(1-fracSens)], 'FontSize', 14)
plot(tlong(end-100), sensfrac_t(end-100), 'k*', 'LineWidth',5)
plot(tlong(end-100), resfrac_t(end-100), 'k*', 'LineWidth',5)
% text(tvec(30), sensfrac_t(30), ['\phi_{sens_t=30hrs}=', num2str(sensfrac_t(30))], 'FontSize', 14)
% text(tvec(30), resfrac_t(30), ['\phi_{res_t=30hrs}=', num2str(resfrac_t(30))], 'FontSize', 14)
% plot(tvec(30), sensfrac_t(30), 'k*', 'LineWidth',5)
% plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
% text(tvec(504), sensfrac_t(504), ['\phi_{sens_t=3WPT}=', num2str(sensfrac_t(504))], 'FontSize', 14)
% text(tvec(504), resfrac_t(504), ['\phi_{res_t=3WPT}=', num2str(resfrac_t(504))], 'FontSize', 14)
% plot(tvec(504), sensfrac_t(504), 'k*', 'LineWidth',5)
% plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
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
dose = 200;
date = '';
[traj2] = comb_data_fxn(traj, filter_criteria, ntreat, dose, date)

%% Write a function to find phi at a certain time
% If repeat dose is within the time frame measured in the first pusle
% treatment
repeat_trts = 1:1:8;
  Cdox = 200;
figure;
for i = 1:length(repeat_trts)
   % Find the Number of sensitive and resisitant cells predicted at start
   % of second treatment
    tdose = repeat_trts(i)*24*7;
    tfirst = 0:1:tdose;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, tfirst, N0, U1, dt, tdrug);
    ind = find(ismember(tfirst, tdose), 1, 'first');
    fracS(i)= Nsri(ind,2)/Nsri(ind,1);
    % Set the initial number of the total cells and the sensitive and
    % resistant cells at the time of repeat treatment
    
    N01(1,i) = N0; % lets say that we simulated reseeding the same number of cells
    % as the first treatment
    % We assume the proportion of sensitive and resistant cells remains the
    % same
    
% Simulate second dose 
    tfin = 672; % this is how long we want to monitor the repeat dose (about 5 weeks)
    tvec2 = 0:dt:tfin;
    U2=k*Cdox*exp(-kdrug*(tvec2))/(0.1*Cdoxmax);
    tdrug = 1; % only simulating one dose)
    % all parameters are the same except now we have a different proportion
    % of sensitive and resistant cells
    p2 = [fracS(i), p(2:end)];
    [Nsr2(:,:,i), tcrit2(i), Ncrit2] = fwd_Greene_model2(p2, tvec2, N0, U2, dt, tdrug);
    Nmodpred = Nsr2(:,:,i);
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
%%
figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,1,i), 'LineWidth', 2)
    hold on
    title('N(t)')
end

figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,2,i), 'LineWidth', 2)
    hold on
    title('S(t)')
end

figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,3,i), 'LineWidth', 2)
    hold on
    title('R(t)')
end
figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,2,i)/Nsr2(:,1,i), 'LineWidth', 2)
    hold on
    title('\phi_{sens}(t)')
end
