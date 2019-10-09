% This script loads in the 231 single pulse treatment parameter estimates
% and makes a prediction for a new dose that wasn't used for calibration.

% We then can compare the prediction of the new dose to the experimental
% data.

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

% Pareto table
pareto_table = readtable('../out/pbest_table.csv');


[phi0, carcapNf, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot] = deal(P{:});

%% Run this forward for a single pulse treatment at 200 nM 
% Again, assume R0=0; dr = 0 and

p = [ phi0, carcapNf,rs,alpha, zr, ds, zd];
% ptest = [pareto_table.phi0vals(m), carcapNf, pareto_table.rsvals(m),...
%     pareto_table.alphavals(m), pareto_table.rr_rs_ratio(m), pareto_table.dsvals(m),...
%     pareto_table.dr_ds_ratio(m)];
pup = [CIphi0(1), carcapNf, CIrs(2), CIalpha(2), CIzr(2), CIds(1), CIzd(1)];
plow = [CIphi0(2), carcapNf, CIrs(1), CIalpha(1), CIzr(1), CIds(2), CIzd(2)];
% Change the CI here to be the min and max values from the table
% other option is to instead sample from the rows of the table, generate
% 100 trajectories, and then just plot the 

% doses used for calibrating
dosevec = [ 1, 4, 7, 9];
Cdoxmax = 1000;
for i= 1:length(trajsum)
Nsri= [];
Nsrlow = [];
Nsrhigh = [];
Nsrmat = [];
dt = 1;
tdrug = 1;
Cdox = trajsum(i).Cdox;
tgen = 0:1:trajsum(i).tvec(end);
Udata=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
N0 = trajsum(i).Nmean(1);
    for m = 1:height(pareto_table)
        pi = [pareto_table.phi0vals(m), carcapNf, pareto_table.rsvals(m),...
            pareto_table.alphavals(m), pareto_table.rr_rs_ratio(m), pareto_table.dsvals(m),...
            pareto_table.dr_ds_ratio(m)];
        [Nsri, tcrit, Ncrit] = fwd_Greene_model2(pi, trajsum(i).tvec, N0, Udata, dt, tdrug);
        Nsrmat = horzcat(Nsrmat,Nsri(:,1));
    end
[Nsrmod, tcrit, Ncrit] = fwd_Greene_model2(p, trajsum(i).tvec, N0, Udata, dt, tdrug);
[Nsrlow, tcrit, Ncritl] = fwd_Greene_model2(plow, trajsum(i).tvec, N0, Udata, dt, tdrug);
[Nsrhigh, tcrit, Ncrith] = fwd_Greene_model2(pup, trajsum(i).tvec, N0, Udata, dt, tdrug);
% Add to the trajsum data structure the model fit for that Cdox (whether it
% is calibrated or predicted)
trajsum(i).Nsrmod = Nsrmod(:,1);
trajsum(i).Nsrmat = Nsrmat;
trajsum(i).Nsrlow = Nsrlow;
trajsum(i).Nsrhigh = Nsrhigh;
end


figure;
for i = 1:4
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    %plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)

    plot(trajsum(j).tvec, trajsum(j).Nsrmod, '-', 'LineWidth', 2,'color', trajsum(j).color )
    for m = 1:height(pareto_table)
    pl2=plot(trajsum(j).tvec, trajsum(j).Nsrmat(:,m), 'k-', 'LineWidth', 2);
    pl2.Color(4) = 0.2;
    end
end
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) data vs model from calibrated doses')
set(gca,'FontSize',20,'LineWidth',1.5)
%%
figure;
for i = 1:4
    subplot(2,2,i)
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    %plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    plot(trajsum(j).tvec, trajsum(j).Nsrmod, 'k-', 'LineWidth', 2 )
    for m = 1:height(pareto_table)
    pl3 = plot(trajsum(j).tvec, trajsum(j).Nsrmat(:,m), 'k-', 'LineWidth', 1);
    pl3.Color(4) = 0.3;
    end
    %plot(trajsum(j).tvec, trajsum(j).Nsrlow(:,1), 'b--', 'LineWidth', 1)
    %plot(trajsum(j).tvec, trajsum(j).Nsrhigh(:,1), 'b--', 'LineWidth', 1)
    xlabel ('time (hours)')
    xlim([0 trajsum(j).tvec(end)])
    
ylabel ('N(t)')
title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM'])
legend('N(t) data', 'model','Location', 'NorthWest')
    legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
end

%%

dosevecnew = [2, 3, 5, 6, 8, 10]
figure;
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)
    plot(trajsum(j).tvec, trajsum(j).Nsrmod(:,1), '-', 'LineWidth', 3, 'color', trajsum(j).color)
    for m = 1:height(pareto_table)
    pl3 = plot(trajsum(j).tvec, trajsum(j).Nsrmat(:,m), 'k-', 'LineWidth', 1);
    pl3.Color(4) = 0.3;
    end
    
    %plot(trajsum(j).tvec, trajsum(j).Nsrlow(:,1), '--', 'LineWidth', 1,'color', trajsum(j).color)
    %plot(trajsum(j).tvec, trajsum(j).Nsrhigh(:,1), '--', 'LineWidth', 1, 'color', trajsum(j).color)
end
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) data vs model predictions at new doses')
set(gca,'FontSize',20,'LineWidth',1.5)

%% subplots of individual doses
Npred = [];
Ndat = [];
figure;
for i = 1:length(dosevecnew)
    subplot(2, 3, i)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    hold on
    plot(trajsum(j).tvec, trajsum(j).Nsrmod(:,1), 'b-', 'LineWidth', 3)
    for m = 1:height(pareto_table)
    pl3 = plot(trajsum(j).tvec, trajsum(j).Nsrmat(:,m), 'k-', 'LineWidth', 1);
    pl3.Color(4) = 0.3;
    end
    
    %plot(trajsum(j).tvec, trajsum(j).Nsrlow(:,1), 'b--', 'LineWidth', 1)
    %plot(trajsum(j).tvec, trajsum(j).Nsrhigh(:,1), 'b--', 'LineWidth', 1)
    % find CCC
    CCC(i) = f_CCC([trajsum(j).Nsrmod(:,1), trajsum(j).Nmean], 0.05);
    Npred = vertcat(Npred, trajsum(j).Nsrmod(:,1));
    Ndat = vertcat(Ndat, trajsum(j).Nmean);
    
    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, CCC= ', num2str(round(CCC(i),3))])
    legend('N(t) data', 'prediction', 'model uncertainty','Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end
CCC_all = f_CCC([Npred, Ndat], 0.05);
figure;
hold on
plot(Ndat, Npred, 'o')
plot([0 max(Npred)], [0 max(Npred)], 'k-', 'LineWidth', 3)
xlim([0 max(Npred)])
ylim([0 max(Npred)])
xlabel('N(t) data')
ylabel('Predicted N(t)')
title(['CCC_{all doses}=', num2str(round(CCC_all,3))])
set(gca,'FontSize',20,'LineWidth',1.5)
%% Now repeat this, but plot all of the raw data for the dose
ct_in_tot = 0;
ct_tot = 0;
figure;
for i = 1:length(dosevecnew)
    subplot(2, 3, i)
    j = dosevecnew(i);
    hold on
    plot(trajsum(j).tvec, trajsum(j).Nsrmod(:,1), 'b-', 'LineWidth', 3)
    plot(trajsum(j).tvec, trajsum(j).Nsrlow(:,1), 'b--', 'LineWidth', 1)
    ct_in = 0;
    for k = 1:6
        plot(trajsum(j).tvec, trajsum(j).Nfit(:,k), '.', 'color', trajsum(j).color)
        ilow=find(trajsum(j).Nfit(:,k)<trajsum(j).Nsrhigh(:,1));
        ihigh = find(trajsum(j).Nfit(:,k)>trajsum(j).Nsrlow(:,1));
        iin=ismember(ihigh,ilow);
        %plot(trajsum(j).tvec(iin), trajsum(j).Nfit(iin,k), 'm*')
        nin = sum(iin);
        ct_in = ct_in+nin;
        ct_totdose= length(trajsum(j).tvec)*6; % total number of data points for that dose
       
    end
    summerr(j) = abs(sum(sum((trajsum(j).Nfit-trajsum(j).Nsrmod(:,1))./trajsum(j).Nfit)));
    pct_in(j) = ct_in/ct_totdose*100;
    MAPE(j) = 100/ct_totdose*summerr(j);
    trajsum(j).pct_in = pct_in(j);
    plot(trajsum(j).tvec, trajsum(j).Nsrhigh(:,1), 'b--', 'LineWidth', 1)
    ct_in_tot = ct_in_tot + ct_in;
    ct_tot = ct_tot+ct_totdose;

    
    % Within this loop, find the number of points that are within the CI
    % bounds and count the percent within
    
    hold on
    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, MAPE= ', num2str(round(MAPE(j),1)), ' %'])
    %title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, ', num2str(round(pct_in(j),0)), '% in bounds'])
    legend('prediction', 'uncertainity',' N(t) data', 'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end
MAPEtotal = 100*sum(summerr)./ct_tot
pct_in_tot = (ct_in_tot/ct_tot)*100
%%

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
