% Run simulations for different dosing regimens using the Greene model
% Start by assuming:
% 1. all cells are sensitive at the start
% 2. resistant cells are invincible (dr = 0)

close all; clear all; clc
%%
dt = 1; % this corresponds to hours
carcap = 2e4; % K
t1 = 336; % 2 weeks
t1vec = (1:dt:t1)';
S = zeros([t1,1]);
R = zeros([t1,1]);
N = zeros([t1,1]);
S(1)=2e3; % set initial conditions (assume all cells sensitive to start)
R(1) = 0; 
Cdox(1) = 75; % nm doxorubicin
rs = 0.01;
rr = 0.001;
ds = 0.0015;
dr = 0;
alpha = 0.0001;
kdrug = 0.025;
N(1) = S(1)+R(1);
t = 1; % time in hours
k = 0.5; 


S0 = S(1);
R0 = R(1);
pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];

% U(t) for a single pulse dose
U1=k*Cdox*exp(-kdrug*(t1vec-1)); % minus one because t starts at 1
tdrug1 =1;

[Nsr, tcrit, Ncrit] = fwd_Greene_model(p, t1vec, U1, dt, tdrug1);
N = Nsr(:,1);
S = Nsr(:,2);
R = Nsr(:,3);
%% Plot for single dose for two weeks of monitoring
figure;
subplot(1,2,1)
plot(t1vec,N,'LineWidth',3, 'color','b');
hold on
plot(t1vec, S, 'LineWidth', 3,'color', 'g')
plot(t1vec, R, 'LineWidth', 3, 'color', 'r')
plot(tcrit, Ncrit, 'k*')
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 t1vec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
title('Single pulse dose response (2 weeks)')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
plot(t1vec, U1,'LineWidth',3)
xlim([ 0 t1vec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
title('Pulse of 75 nM dox')
%% Simulate second dose at two weeks

t2vec = (1:1:336)'; %
tvec = vertcat(t1vec, t2vec + t1vec(end));
tdrug = [t1vec(1); (t2vec(1) + t1vec(end))];
% set original U as U1 to indicate first dose

U2 = k*Cdox*exp(-kdrug*(t2vec-1)); % minus one because t starts at 1
U = vertcat(U1, U2);
[Nsr, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
N = Nsr(:,1);
S = Nsr(:,2);
R = Nsr(:,3);

figure;
subplot(1,2,1)
plot(tvec,N,'LineWidth',3, 'color','b');
hold on
plot(tvec, S, 'LineWidth', 3,'color', 'g')
plot(tvec, R, 'LineWidth', 3, 'color', 'r')
for i = 1:length(tdrug)
plot(tcrit(i) + tdrug(i)-1, Ncrit(i), 'k*', 'LineWidth',3)
text(tcrit(i)+tdrug(i)+2, Ncrit(i), ['t_{crit}=', num2str(tcrit(i)), ' hrs post treat'])
end
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
title('Multiple treatment response')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
plot(tvec, U,'LineWidth',3)
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
title('Dosing regimen')

%% Simulate multiple doses
% Set up easy way to make time and dose vector
int_treat = [2; 2; 2]; % intervals between treatments
cum_treat = cumsum(int_treat);
totweeks = sum(int_treat) + 4; % monitor for two weeks after last treatment
tdrug = [1; cum_treat*24*7];
tvec = [1:1: totweeks*7*24]';
Cdox = [ 75; 75; 75; 75];
acc_dose = cumsum(Cdox)
Uvec = [];
Utot = zeros([length(tvec),1]);
for i = 1:length(tdrug)
% start with zeros 1:tdrug(1)
Uinit = zeros([tdrug(i)-1,1]);
tin = tvec(tdrug(i):end)-tdrug(i);
Udrug = k*Cdox(i)*exp(-kdrug*(tin-1)); % minus one because t starts at 1
% add on each additional treatment
Upulse = vertcat(Uinit, Udrug);
Utot = Utot + Upulse;

end
alpha = 1e-4;
pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];


[Nsr, tcrit, Ncrit] = fwd_Greene_model(p, tvec, Utot, dt, tdrug);
N = Nsr(:,1);
S = Nsr(:,2);
R = Nsr(:,3);

figure;
subplot(1,3,1)
plot(tvec,N,'LineWidth',3, 'color','b');
hold on
plot(tvec, S, 'LineWidth', 3,'color', 'g')
plot(tvec, R, 'LineWidth', 3, 'color', 'r')
for i = 1:length(tdrug)
plot(tcrit(i) + tdrug(i)-1, Ncrit(i), 'k*', 'LineWidth',3)
text(tcrit(i)+tdrug(i)+2, Ncrit(i), ['t_{crit}=', num2str(tcrit(i)), ' hrs post treat'])
end
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
title('Multiple treatment response')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,3,2)
plot(tvec, Utot,'LineWidth',3)
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose (U(t))','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
title('Dosing regimen')

subplot(1,3,3)
plot(acc_dose, tcrit, 'k*', 'LineWidth',3)
xlabel('Accumulated Dose (nM)','FontSize',20)
ylabel('T_{crit} post last treatment','FontSize',20)
title(['T_{crit} vs. Accumulated Dose. alpha= ',num2str(alpha)])
set(gca,'FontSize',20,'LineWidth',1.5)