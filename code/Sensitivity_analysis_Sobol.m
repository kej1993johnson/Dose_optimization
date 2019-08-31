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
S = load('../out/trajfit231.mat');
traj= S.traj;
Ssum = load('../out/trajsumfit231.mat');
trajsum = Ssum.trajsum;

% Use the fit parameters from N(t) only to make the parameter domains
% We shouldn't really need these... but maybe keep just as a test
ptest = load('../out/ptest.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);
% Can use some of these as first guesses/ballpark ranges of what to expect
[phi0f, carcapNf, carcapphif, rsf, alphaf, zrdata, dsf, zdf, k, kdrug, gtot]= deal(P{:});
%[rsf, carcapNf, alphaf, rrf, dsf, drf, k, kdrug, gtot] = deal(P{:});
% We will use this carcapNf only, and the k and kdrug to be consistent

phi_est_filename = '../data/phi_t_est.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t;
phitrt = phi_est.phi_t;
tscRNAseq = tbot(end);

%% Plot the average data for the set of controls (u(t)s) we want to use
  figure;
 for i = 1:length(trajsum)-4
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

%% Run loop to generate output needed from each control u(t)
% If we want to output the effect of individual parameter changes on the
% "model behavior" and we quantify model behavior as some observable
% feature about the treatment response-- then we could say let's have this
% output tcrit versus dox concentration & phi(tseq) versus treatment (i.e.
% phenotypic composition at the time of treatment
Cdoxmax = 1000;
Cdoxvec = [0, 25, 50, 75];
Cdoxvec =[0:10:150];
tvec1=[0:3:450];
tgen1 =[0:1:tvec1(end)];
tvec2 =[0:3:tscRNAseq];
tgen2 = [0:1:tscRNAseq];
tdrug = 1;
for i = 1:length(Cdoxvec)
U1(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen1))/(0.1*Cdoxmax);
dt = 1;
tdrug = 1;
N0 = 2e3;
p1 = [ptest(1:2), ptest(4:8)];
% Use this to find the critical times as a function of dose in 96 well
[Nsrdat1(:,:,i), tcrit1(i), Ncrit1(i)]=fwd_Greene_model2(p1, tvec1, N0, U1(:,i), dt,tdrug);
% Use this to find the phi(@tscRNAseq) as a function of dose in expansion
tgen2 = [0:1:tscRNAseq];
U2(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen2))/(0.1*Cdoxmax);
p2 =p1;
% replace carcaph
p2(2) = carcapphif;
[Nsrdat2(:,:,i), tcrit2(i), Ncrit2(i)]=fwd_Greene_model2(p2, tvec2, N0, U2(:,i), dt,tdrug);
% Use this to find phi@tSCRNAseq
% This could be found at any time, we just pick last bc last is the time of
% scNAseq
phi_sc(i) = Nsrdat2(end,2,i)/Nsrdat2(end,1,i);
end
%% Plot an example and then the resulting tcrit and phi(scRNAseq) vs. Cdox.
figure;
subplot(1,3,1)
plot(tgen1, U1(:,2), 'k', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Effective dose (u(t))')
title(['Dosing input, Cdox= ', num2str(Cdoxvec(2)), ' nM'])
xlim([0, tgen1(end)])

subplot(1,3,2)
plot(tvec1, Nsrdat1(:,1,2), 'b', 'LineWidth', 2)
hold on
plot(tvec1,Nsrdat1(:,2,2), 'g', 'LineWidth', 2)
plot(tvec1, Nsrdat1(:,3,2), 'r', 'LineWidth',2)
plot(tcrit1(2), Ncrit1(2), 'k*', 'LineWidth', 2)
text(tcrit1(2), Ncrit1(2), ['t_{crit}= ',num2str(tcrit1(2)),' hours'], 'FontSize',14,'HorizontalAlignment', 'left')
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('N(t)')
title('N(t) from 96 well')
xlim([0, tgen1(end)])

subplot(1,3,3)
plot(tvec2, Nsrdat2(:,1,2), 'b', 'LineWidth', 2)
hold on
plot(tvec2,Nsrdat2(:,2,2), 'g', 'LineWidth', 2)
plot(tvec2, Nsrdat2(:,3,2), 'r', 'LineWidth',2)
text(tscRNAseq, Nsrdat2(end,1,2), ['\phi(@scRNAseq)=',num2str(phi_sc(2))], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Number of cells')
title('N(t) from scRNAseq')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, tgen2(end)])

figure;
subplot(1,2,1)
plot(Cdoxvec, tcrit1, 'b-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('Dox concentration')
ylabel('Time to reach 2N_{0}')
title('Dox vs. t_{crit}')
xlim([0 Cdoxvec(end)])
subplot(1,2,2)
plot(Cdoxvec, phi_sc, 'm-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('Dox concentration')
ylabel('\phi_{s} at time of scRNAseq (10 weeks)')
title('Dox vs. \phi_{s} @ 10 weeks')
xlim([0 Cdoxvec(end)])
%% Perform a local sensitivity analysis
% tornado_plot of sensitivities of model at the default values for other parameters
 

par = ptest(1:8);
parlow = par/2;
parhigh = par*2;

paramnames = {'\phi_{s}(0)', 'K1', 'K2', 'r_{s}', '\alpha', '(rr/rs)', 'd_{s}', '(dr/ds)'};

Npar = length(par);

% Write in line functions that output tcrit and phi_sc 
modelfuntcrit=@(p)drug_ind_wrapper(p,Cdoxvec);
modelfunphisc=@(p)drug_ind_wrapper2(p,Cdoxvec);

defvaltcrit = modelfuntcrit(par)';
defvalphisc = modelfunphisc(par)';

for j = 1:Npar
	plow = par; phigh = par;
    % Replace each parameter one by one with high and low value
	plow(j)=parlow(j);
    % record the result of lowering the jth parameter
	lovaltcrit(:,j) = modelfuntcrit(plow)';
    lovalphisc(:,j) = modelfunphisc(plow)';
    % Replace jth parameter with high value
	phigh(j) = parhigh(j);
    % record the result of raising the jth parameter
	hivaltcrit(:,j) = modelfuntcrit(phigh)';
    hivalphisc(:,j) = modelfunphisc(phigh)';
end

%% now plot the results for tcrit
defval = defvaltcrit;
loval = lovaltcrit;
hival = hivaltcrit;

figure;
sq_diffs_low = (defval-loval).^2;
sq_diffs_high = (hival - defval).^2;
diffs_lows = sum(sq_diffs_low,1);
diffs_high = sum(sq_diffs_high,1);
diffs = diffs_high + diffs_lows ; 
[~,jsort] = sort(diffs);
def_vals_zeros = zeros(1,Npar);
for j = 1:Npar
	%plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot([-diffs_lows(jsort(j)) diffs_high(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
    plot(-diffs_lows(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(diffs_high(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(def_vals_zeros,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if diffs_lows(jsort(j)) < diffs_high(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	%text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	%text(diffs_lows(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'low') '  '],'HorizontalAlignment',loalign,'color','g', 'FontSize',12);
	%text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'high') '  '],'HorizontalAlignment',highalign,'color','r', 'FontSize',12);
    %text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j)))) '  '],'HorizontalAlignment','right','color','k', 'FontSize',12);
end

hold on
plot(0*[1 1],[0 Npar+0.5],'k-');
%set(gca,'Ytick',[1:Npar],'Yticklabel',par(jsort), 'FontSize', 15);
set(gca,'Ytick',[1:Npar],'Yticklabel',char(paramnames(jsort)), 'FontSize', 15);
%set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('t_{crit}');
grid on
title('Local Sensitivity in t_{crit}')
set(gca,'FontSize',16)
%% Now for phisc
defval = defvalphisc;
loval = lovalphisc;
hival = hivalphisc;
figure;
sq_diffs_low = (defval-loval).^2;
sq_diffs_high = (hival - defval).^2;
diffs_lows = sum(sq_diffs_low,1);
diffs_high = sum(sq_diffs_high,1);
diffs = diffs_high + diffs_lows ; 
[~,jsort] = sort(diffs);
def_vals_zeros = zeros(1,Npar);
for j = 1:Npar
	%plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot([-diffs_lows(jsort(j)) diffs_high(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
    plot(-diffs_lows(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(diffs_high(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(def_vals_zeros,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if diffs_lows(jsort(j)) < diffs_high(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	%text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	%text(diffs_lows(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'low') '  '],'HorizontalAlignment',loalign,'color','g', 'FontSize',12);
	%text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'high') '  '],'HorizontalAlignment',highalign,'color','r', 'FontSize',12);
    %text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j)))) '  '],'HorizontalAlignment','right','color','k', 'FontSize',12);
end

hold on
plot(0*[1 1],[0 Npar+0.5],'k-');
%set(gca,'Ytick',[1:Npar],'Yticklabel',par(jsort), 'FontSize', 15);
set(gca,'Ytick',[1:Npar],'Yticklabel',char(paramnames(jsort)), 'FontSize', 15);
%set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('\phi_{s} @ 10 weeks');
grid on
title('Local Sensitivity in \phi_{s}')
set(gca,'FontSize',16)
%% Next step is to use the modelfuns to perform a global sensitivity analysis
% Using the procedure outlined in Angela's slides

