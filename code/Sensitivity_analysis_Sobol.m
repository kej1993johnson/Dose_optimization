% This script is meant to be used to perform a global sensitivity analysis
% on the parameter ranges derived from the calibrated parameters.  The
% results of this sensitivity analysis will reveal the most important
% parameters of the system, causing the greatest variation in outputs, for
% the area of parameter space for which the model is able to replicate the
% experimental data and the area of uncertainty.

% The important parameters have been shown to be drivers in changing the
% stability of steady states in mathematical models and may drive future
% experimental investigations and/or provide support to theories about the
% biology. Could also be used to reveal potential simplified versions of
% the model. 

% We start by doing this on the simplest model of:
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R

% For now start by keeping the u(t) parameters out of the analysis, but can
% put these in if we need to later

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/trajfit.mat');
traj= S.traj;

S = load('../out/trajsumfit.mat');
trajsum = S.trajsum;

p4fit = load('../out/p4fit.mat');
p4fit = struct2cell(p4fit);
p4fit = cell2mat(p4fit);
P = num2cell(p4fit);
[rs, carcap, alpha, rr, ds, dr] = deal(P{:});
%% Run this forward for a single pulse treatment at 75 nM 
% Again, assume R0=0; dr = 0 and
kdrug = 0.0175;
k = 0.5;
tdrug =1; % stands for first treatment
dt = 1;


% Pull from trajsum dose = 75
for i = 1:length(trajsum)
    if trajsum(i).Cdox == 75
        tvec = trajsum(i).tvec;
        Nmean = trajsum(i).Nmean;
        Nstd = trajsum(i).Nstd;
    end
end
tvec = 0:1:504; % three weeks 
R0 = 0;
S0 = Nmean(1);
p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];

Cdox = 75;
U1 = k*Cdox*exp(-kdrug*(tvec));
[Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U1, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
Nsrend = Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);

figure;
plot(tvec, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tvec, Nsri(:,2), 'g','LineWidth',2)
plot(tvec, Nsri(:,3),'r', 'LineWidth',2)
text(tvec(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)])
text(tvec(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)])
text(tvec(end-20), Nsri(end-20,1),['N{_f}=', num2str(Nsrend)])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 75 nM pulse treatment')
legend ('N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
%% Run loop to generate output needed
p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
N0 = Nmean(1);
% First need to establish a domain for every parameter


% First do a local sensitivity analysis
[Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U1, dt, tdrug);
Nsrendcurr = Nsri(end,1);
fracSensendcurr = Nsri(end,2)./Nsri(end,1);
Nsrend = [];
for i = 1%1:length(p)-1
    qvec = [1,3:length(p)];
    pit = p;
        if i ==1
            pit(1) = p(1)*0.95;
            pit(2) = N0-pit(1);
            deltap = pit(1)-p(1);
        elseif i == 7
            pit(7) = p(7) + 1e-10;
            deltap = 1e-10;
        elseif i~=1
        pit(qvec(i))=p(qvec(i))*0.95;
        deltap = pit(qvec(i))-p(qvec(i));
        end
        Cdox = 75;
        U1 = k*Cdox*exp(-kdrug*(tvec));
[Nsri, tcrit, Ncrit] = fwd_Greene_model(pit, tvec, U1, dt, tdrug);
fracSensend(i) = Nsri(end,2)./Nsri(end,1);
Nsrend(i) = Nsri(end,1);
sens(1,i) = abs((Nsrend(i)-Nsrendcurr)./(deltap));
sens(2,i) = abs((fracSensend(i)-fracSensendcurr)./(deltap));
end

params = ({'S{0}', 'r_{s}', 'K', '\alpha', 'r_{r}', 'd_{s}', 'd_{r}'});

figure;
bar(1:1:length(p)-1, sens(1,:))
xticks(1:length(p)-1)
xticklabels(params)
ylabel('Local Sensitivity to N_{t=tfinal}')
xlabel('parameters')
title('Local Parameter Sensitivity on N')

figure;
bar(1:1:length(p)-1,sens(2,:))
xticks(1:length(p)-1)
xticklabels(params)
ylabel('Local Sensitivity to \phi_{sens}')
xlabel('parameters')
title('Local Parameter Sensitivity on \phi_{sens}')

%% Find total variance in N and phi
n=10;
S0dom = linspace(1,N0,n);
rsdom = linspace(0,0.01, n);
carcapdom= linspace(2e3,1e6, n);
alphadom = linspace(0, 0.05, n);
rrdom = linspace(0, 0.005, n);
dsdom = linspace(0,0.01,n);
drdom = linspace(0,0.001,n);
mat = vertcat(S0dom, rsdom, carcapdom, alphadom, rrdom, dsdom, drdom);



Nall = reshape(Nsrend, [1,n*(length(p)-1)]);
fracall = reshape(fracSensend, [1, n*(length(p)-1)]);

VarN = var(Nall);
Varphi = var(fracall);
for i = 1:length(p)-1
    Si(1,i) = (var(Nsrend(i,:)))./(VarN);
    Si(2,i) = (var(fracSensend(i,:)))./(Varphi);
end
params = {'S{0}', 'r_{s}', 'K', '\alpha', 'r_{r}', 'd_{s}', 'd_{r}'}
figure;
bar(Si', params)