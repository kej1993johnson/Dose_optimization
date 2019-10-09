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

% Also don't really need the confidence intervals
CI = load('../out/CIpbest.mat');
CI = struct2cell(CI);
CI = cell2mat(CI);
CIphi0 = CI(1,:);
CIrs = CI(2,:);
CIalpha = CI(3,:);
CIzr = CI(4,:);
CIds = CI(5,:);
CIzd=CI(6,:);

[phi0, carcapNf, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot] = deal(P{:});

%% Run this forward for a single pulse treatment at 200 nM 
% Again, assume R0=0; dr = 0 and


% Pull from trajsum dose = 200
for i = 1:length(trajsum)
    if trajsum(i).Cdox == 200
        tvec = trajsum(i).tvec;
        Nmean = trajsum(i).Nmean;
        Nstd = trajsum(i).Nstd;
        N0 = Nmean(1);
    end
end

p = [ phi0, carcapNf,rs,alpha, zr, ds, zd];
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
title('N(t) for single 200 nM pulse treatment')
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
xlim([ 0 384])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 200 nM pulse treatment')
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
xlim([ 0 384])
ylabel ('Proportion of cells')
title('\phi_{S} & \phi_{R} following single 200 nM pulse treatment')
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
repeat_trts = [ 4 6 8 10 ]; % days at which the treatment was repeated
  Cdox = 200;
  Cdoxmax = 1000;
figure;
for i = 1:length(repeat_trts)
   % Find the Number of sensitive and resisitant cells predicted at start
   % of second treatment
    tdose = repeat_trts(i)*24;
    tfirst = 0:1:tdose;
    tsecond = 0:1:336; % monitor for two weeks after
    tdrug = [0; tdose];
    ttot = horzcat(tfirst, tsecond(2:end) + tfirst(end));
    U1=k*Cdox*exp(-kdrug*(tfirst))/(0.1*Cdoxmax);
    U2=k*Cdox*exp(-kdrug*(tsecond))/(0.1*Cdoxmax);
    U = horzcat(U1(1:end-1), U2);
    % Run the model forward for the two different U(t)s and times of drug
    % for each variest treatment interval
    [Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, ttot, N0, U, dt, tdrug);
    

subplot(2, length(repeat_trts)/2, i)

    hold on
    plot (ttot, Nsri(:,1),'b-' ,'LineWidth',2)
    plot(ttot, Nsri(:,2), 'g-', 'LineWidth',2)
    plot(ttot, Nsri(:,3), 'r-', 'LineWidth', 2)
    plot(tvec, Nmean,'k.', 'LineWidth', 2)
    %plot(tvec, Nmean + Nstd, 'color', 'k')
    %plot(tvec, Nmean - Nstd, 'color', 'k')
    title (['Treatment Interval=', num2str(repeat_trts(i)), ' days'])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    legend ('N(t)', 'S(t)', 'R(t)','N(t) data', '95% CI', 'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    xlim([0 ttot(end)])
    ylim([0 0.7*carcapNf])
end

%% Repeat this but add in a cone of uncertainty around N...

% Set the upper bound parameters
%p = [ phi0, carcapNf,rs,alpha, zrdata, ds, zd];
sigtech = 1e-2;
% Combine all the parameters that should lead to a lower response...
pup = [CIphi0(1), carcapNf, CIrs(2), CIalpha(2), CIzr(2), CIds(1), CIzd(1)];
plow = [CIphi0(2), carcapNf, CIrs(1), CIalpha(1), CIzr(1), CIds(2), CIzd(2)];
%

repeat_trts = [4 6 8 10]; % days at which the treatment was repeated
  Cdox = 200;
  Cdoxmax = 1000;
figure;
for i = 1:length(repeat_trts)
   % Find the Number of sensitive and resisitant cells predicted at start
   % of second treatment
    tdose = repeat_trts(i)*24;
    tfirst = 0:1:tdose;
    tsecond = 0:1:336; % monitor for two weeks after
    tdrug = [0; tdose];
    ttot = horzcat(tfirst, tsecond(2:end) + tfirst(end));
    U1=k*Cdox*exp(-kdrug*(tfirst))/(0.1*Cdoxmax);
    U2=k*Cdox*exp(-kdrug*(tsecond))/(0.1*Cdoxmax);
    U = horzcat(U1(1:end-1), U2);
    % Run the model forward for the two different U(t)s and times of drug
    % for each variest treatment interval
    [Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, ttot, N0, U, dt, tdrug);
    % Now simulate upper and lower bounds from parameter estimates
    [Nsrlow, tcritl, Ncritl] = fwd_Greene_model2(plow, ttot, N0, U, dt, tdrug);
    [Nsrhigh, tcrith, Ncrith] = fwd_Greene_model2(pup, ttot, N0, U, dt, tdrug);
subplot(2, length(repeat_trts)/2, i)

    hold on
    plot(ttot, Nsri(:,1),'b-' ,'LineWidth',2)
    plot(ttot, Nsri(:,2),'g-' ,'LineWidth',2)
    plot(ttot, Nsri(:,3),'r-' ,'LineWidth',2)
    
    plot(ttot, Nsrlow(:,1), 'b--', 'LineWidth',1)
    plot(ttot, Nsrhigh(:,1), 'b--', 'LineWidth', 1)
    
%     plot(ttot, Nsrlow(:,2), 'g--', 'LineWidth',1)
%     plot(ttot, Nsrhigh(:,2), 'g--', 'LineWidth',1)
%     
%     plot(ttot, Nsrlow(:,3), 'r--', 'LineWidth',1)
%     plot(ttot, Nsrhigh(:,3), 'r--', 'LineWidth',1)
%     
    
    title (['Treatment Interval=', num2str(repeat_trts(i)), ' days'])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    legend ('N(t)', 'lower N(t)', 'upper N(t)','N(t) data', '95% CI', 'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    xlim([0 ttot(end)])
    %ylim([0 carcapNf])
    ylim([0 0.7*carcapNf])
end
figure;
plot(ttot, U, '-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time(hours)')
ylabel('U(t)')
title('Repeat doses: U(t)')

%% Just add uncertainty to the second treatment.
% Need to figure out why this isnot continuous even though the one calling
% the function is!
figure;
for i = 1:length(repeat_trts)
    
    % Simulate the first dose with the best fitting p only
    tdose = repeat_trts(i)*24;
    tfirst = 0:1:tdose;
    tsecond = 0:1:336; % monitor for two weeks after
    tdrug = [0; tdose];
    ttot = horzcat(tfirst, tsecond(2:end) + tfirst(end));
    U1=k*Cdox*exp(-kdrug*(tfirst))/(0.1*Cdoxmax);
    U2=k*Cdox*exp(-kdrug*(tsecond))/(0.1*Cdoxmax);
    U = horzcat(U1(1:end-1), U2);
    
    % Run the model forward for the two different U(t)s and times of drug
    % for each variest treatment interval
    [Nsr1, tcrit, Ncrit] = fwd_Greene_model2(p, tfirst, N0, U1(1:end-1), dt, tdrug(1));
    
    
    % Now simulate the secodn treatment with
    % upper and lower bounds from parameter estimates
    p2 = p;
    p2(1) = Nsr1(end,2)./Nsr1(end,1);
    plow2 = plow;
    plow2(1) = p2(1) + sigtech;
    pup2 = pup;
    pup2(1) = p2(1)-sigtech;
    [Nsr2, tcrit, Ncrit] = fwd_Greene_model2(p2, tsecond, Nsr1(end,1), U2, dt, tdrug(1));
    [Nsr2low, tcritl, Ncritl] = fwd_Greene_model2(plow2, tsecond, Nsr1(end,1), U2, dt, tdrug(1));
    [Nsr2high, tcrith, Ncrith] = fwd_Greene_model2(pup2, tsecond, Nsr1(end,1), U2, dt, tdrug(1));
   % Ntot = vertcat(Nsr1(1:end,:), Nsr2(:,:));
    iend = find(tvec>tfirst(end), 1, 'first');
    
    subplot(2, length(repeat_trts)/2, i)

    hold on
    plot(tfirst, Nsr1(:,1), 'b-', 'LineWidth', 2)
    plot(tfirst, Nsr1(:,2), 'g-', 'LineWidth', 2)
    plot(tfirst, Nsr1(:,3), 'r-', 'LineWidth',2)
    plot(tvec(1:iend), Nmean(1:iend),'k.', 'LineWidth', 2)
    plot(tvec(1:iend), Nmean(1:iend) + Nstd(1:iend), 'color', 'k')
    plot(tvec(1:iend), Nmean(1:iend) - Nstd(1:iend), 'color', 'k')
    plot(tsecond + tfirst(end), Nsr2(:,1), 'b-', 'LineWidth',2)
    plot(tsecond + tfirst(end), Nsr2(:,2), 'g-', 'LineWidth',2)
    plot(tsecond + tfirst(end), Nsr2(:,3), 'r-', 'LineWidth',2)
    plot(tsecond+tfirst(end), Nsr2low(:,1), 'b--', 'LineWidth',1)
    plot(tsecond+tfirst(end), Nsr2high(:,1), 'b--', 'LineWidth', 1)
    title (['Treatment Interval=', num2str(repeat_trts(i)), ' days'])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    legend ('N(t)', 'S(t)', 'R(t)', 'N(t) data',  'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    xlim([0 ttot(end)])
    %ylim([0 carcapNf])
    ylim([0 0.7*carcapNf])
    
end


%%
figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,1,i), 'LineWidth', 2)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('N(t) for different repeat treatment intervals' )
    xlim([0 tfin])
end
set(gca,'FontSize',20,'LineWidth',1.5)
legend('1 WPT', '2 WPT', '3 WPT', '4 WPT', '5 WPT', '6 WPT', '7 WPT','8 WPT', 'Location', 'Northwest')
legend boxoff

figure;
for i = 1:8
    plot(tvecmod, Nsr2(:,2,i)/Nsr2(:,1,i), 'LineWidth', 2)
    hold on
    xlabel('time (hours)')
    ylabel('\phi_{s}(t)')
    title('\phi_{s}(t) for different repeat treatment intervals' )
     xlim([0 tfin])
end
set(gca,'FontSize',20,'LineWidth',1.5)
legend('1 WPT', '2 WPT', '3 WPT', '4 WPT', '5 WPT', '6 WPT', '7 WPT','8 WPT', 'Location', 'Northwest')
legend boxoff

%%
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
