% Make dose history figure

% This script will be used to make visualizations of the dosing history for
% each dose response curve in the data set, using the information in traj.
% mat.

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajraw.mat');
traj= S.traj;
%% Start with an example for 2 repeat treatments 3 weeks apart at 75 nM each
% this corresponds to the 32th row in the structure
i = 32;

% Start my making a time vector that spans the pre-treatment history + two
% weeks into the current treatment. So here that will be
tvec = 0:1:(7*(sum(traj(i).doseints)+2));

% Next make your dose vec. Start with a vector of zeros the length of your
% time vector
Uvec = 0.*tvec;
% Add in pulse treatments at the corresponding days
% First pulse treatment starts at 0
ind1 = find(ismember(tvec, 0:traj(i).doseduration/24))
Uvec(ind1) = traj(i).prevdose(1)
%% Second pulse treatment start at variable places given by the 

for j = 1:length(traj(i).doseints)
    ind = find(ismember(tvec, (traj(i).doseints(j)*7):((traj(i).doseints(j)*7)+traj(i).doseduration/24)));
    Uvec(ind)= traj(i).prevdose(j);
end

figure;
plot(tvec, Uvec, 'b-', 'LineWidth', 3)
xlabel('time (days)')
ylabel('Dose')
title('Treatment History Example')
set(gca,'FontSize',20,'LineWidth',1.5)
