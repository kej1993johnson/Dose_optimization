% Make dose history figure

% This script will be used to make visualizations of the dosing history for
% each dose response curve in the data set, using the information in traj.
% mat.

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajraw.mat');
traj= S.traj;
%% pick well (row in traj) to graph
i = 58;


%% following block is Kaitlyn's old time and dosage vectors
%{
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
%}
%% time and dosage vectors to plot
% generates tvec truncated after last dose
tvec = 0:1:(sum(traj(i).doseintdays) + 2 + (3 * traj(i).numdoses));
% equalize vector length of y values
uvec = 0.*tvec;
n = traj(i).numdoses;
% see "dosegraphindex.pdf" in documentation to understand this index jumping/modifying
ind = 4;
for j = 1:n
    for k = ind:length(tvec)
        tvec(k) = tvec(k) - 1;
    end
    for k = (ind + 2):length(tvec)
        tvec(k) = tvec(k) - 1;
    end
    uvec(ind) = traj(i).dose;
    uvec(ind + 1) = traj(i).dose;
    if j ~= n
        ind = ind + 3 + traj(i).doseintdays(j);
    end
end    
%% plot and format
plot(tvec, uvec)
uvec;
xlabel('time (days)')
ylabel('Dose')
title('Treatment History Example')
set(gca,'FontSize',20,'LineWidth',1.5)
