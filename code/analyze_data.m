% This script analyzes the MCF-7 data-- lot that can be done here
% So far just wanted to plot some things

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajfit.mat');
traj= S.traj;
%% Plot some things that might be interesting
% Just playing around here

figure;
subplot(1,3,1)
for i = 1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).phi, 'r*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('resistant fraction (\phi)')
    title('Resistant fraction versus WPT')
    end
end
subplot(1,3,2)
for i =1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).k, 'b*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('death rate (k)')
    title('Death rate versus WPT')
    end
end

subplot(1,3,3)
for i =1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).g, 'g*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('resistant regrowth rate (g)')
    title('Resistant regrowth rate versus WPT')
    end
end