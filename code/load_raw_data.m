% This script loads in data sets from Grant's various experimental
% conditions of drug response. The goal of this script is just to capture
% the input variables and data for each well in the data set.

close all; clear all; clc

% load in adjusted data set (remove top part and record relevant info
% below)
[raw(1).N, raw(1).T] =xlsread('../data/GH1831_MCF7_PredModel_2nd_Dose_Varying_Interval_v2.xls');
[raw(2).N, raw(2).T] =xlsread('../data/GH1830_MCF7_PredModel_3rd_Dose_Vary_Interval_v2.xls');
[raw(3).N, raw(3).T] =xlsread('../data/GH1818_MCF7_PredModel_1st_Dose_Analysis_v6.xls');
[raw(4).N, raw(4).T] =xlsread('../data/GH1832_MCF7_PredModel_3rd_Dose_4-4_Intervals_v1.xls');
[raw(5).N, raw(5).T] =xlsread('../data/GH1820_MCF7_PredModel_2nd_Dose_Varying_Interval_v3.xls');
[raw(6).N, raw(6).T] =xlsread('../data/GH1829_MCF7_PredModel_2nd_Dose_4_Week_Interval_v1.xls');
[raw(7).N, raw(7).T] =xlsread('../data/GH1819_MCF7_PredModel_3rd_Treatment_Analysis_v2.xls');
[raw(8).N, raw(8).T] =xlsread('../data/GH1817_MCF7_PredModel_2nd_Dose_Analysis_v3.xls');
[raw(9).N, raw(9).T] =xlsread('../data/GH1815_MCF7_PredModel_Effect_of_1st_Dose_Analysis_v2.xls');
[raw(10).N, raw(10).T] =xlsread('../data/GH1835_BT474_PredModel_2nd_Dose_Vary_Interval_v1.xls');

% this reads in the numeric values in an excel sheet (N1), and the text
% values (T1)
% Make a "size" matrix that keeps track of the number of wells and time
% points of each data set loaded. This will be used to add to the structure
% once more data sets are loaded in
% 1st column in number of time points, 2nd column number of wells in that
% data set, third column is cumulative number of wells
<<<<<<< HEAD
sz(1,1:2) = size(N1);
sz(1,2)=sz(1,2)-1;

sz(2,1:2) = size(N2);
sz(2,2) = sz(2,2)-1;

sz(3,1:2) = size(N3);
sz(3,2) = sz(3,2)-1;

sz(4,1:2) = size(N4);
sz(4,2) = sz (4,2) - 1;

sz(5,1:2) = size(N5);
sz(5,2) = sz (5,2) - 1;

sz(6,1:2) = size(N6);
sz(6,2) = sz (6,2) - 1;

sz(7,1:2) = size(N7);
sz(7,2) = sz (7,2) - 1;

sz(8,1:2) = size(N8);
sz(8,2) = sz(8,2) - 1;

sz(9,1:2) = size(N9);
sz(9,2) = sz (9,2) - 1;

sz(:,3) = cumsum(sz(:,2));

raw(1).N = N1;
raw(2).N = N2;
%% Make structuren (1831)
%Make one structure where each entry is a trajectory containing:
% 1. time vector
% 2. cell  number in time
% 3. date of dosing
% 4. dose doxorubicin (nM)
% 5. dose number (i.e. 1st, 2nd, 3rd)
% 6. time after previous dose (in weeks)
% 7. previous doses (nM)
% 8. intervals between doses
% 9. accumulated dose
% 10. number of doses
% 11. seeding density (cells/well)

for i = 1:sz(1,2) % size matrix, first row, second column
    traj(i).time = N1(1:end,1);
    traj(i).rawN = N1(1:end,i+1);
    traj(i).date = '12-19-18'; % get this from original excel file
    traj(i).welllabel = T1(1, i+1);
    welllabel = string(T1(1, i+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose =  60;
 
    
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
    cond8 = { 'B2', 'C2', 'D2','E2', 'F2','G2'};
    cond7 = { 'B3', 'C3', 'D3','E3', 'F3','G3'};
    cond6 = { 'B4', 'C4', 'D4','E4', 'F4','G4'};
    cond5 = { 'B5', 'C5', 'D5','E5', 'F5','G5'};
    cond4 = { 'B6', 'C6', 'D6','E6', 'F6','G6'};
    cond3 = { 'B7', 'C7', 'D7','E7', 'F7','G7'};
    cond2 = { 'B8', 'C8', 'D8','E8', 'F8','G8'};
    cond1 = { 'B9', 'C9', 'D9','E9', 'F9','G9'};
    cond5XT = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    condUT = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
    if contains(traj(i).well, cond8 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 8;
        traj(i).prevdose = 75; % note if more than two treatments this could be a vector
        traj(i).doseints = 8; % same with this 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000; % this is from the number of cells intended to be seeded (from well label)
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond7 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 7;
        traj(i).prevdose = 75; 
        traj(i).doseints = 7;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond6 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 6;
        traj(i).prevdose = 75; 
        traj(i).doseints = 6;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond5 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 5;
        traj(i).prevdose = 75; 
        traj(i).doseints = 5;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond4 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 4;
        traj(i).prevdose = 75; 
        traj(i).doseints = 4;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond3 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 3;
        traj(i).prevdose = 75; 
        traj(i).doseints = 3;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond2 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 2;
        traj(i).prevdose = 75; 
        traj(i).doseints = 2;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond1 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 1;
        traj(i).prevdose = 75; 
        traj(i).doseints = 6;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end

    if contains(traj(i).well, cond5XT )
        traj(i).dose = 75;
        traj(i).dosenum = 6;
        traj(i).WPT = 8;
        traj(i).prevdose = [75, 100, 75, 75, 75]; 
        traj(i).doseints = [2, 4, 3, 4, 8];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 6;
        traj(i).seed = 2000;
        traj(i).flag = 1;
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 75;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
end
%% Load in second data set (1830: 3rd dose with different intervals)

for i = sz(1,3)+1:sz(2,3) % size matrix, first row, second column
    k=i-sz(1,3);
    traj(i).time = N2(1:end,1);
    traj(i).rawN = N2(1:end,k+1);
    traj(i).date = '12-19-18'; % get this from original excel file
    traj(i).welllabel = T2(1, k+1);
    welllabel = string(T2(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 69;
 
    
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
    cond4_4 = { 'B2', 'C2', 'D2','E2', 'F2','G2'};
    cond4_3 = { 'B3', 'C3', 'D3','E3', 'F3','G3'};
    cond3_4 = { 'B4', 'C4', 'D4','E4', 'F4','G4'};
    cond4_2 = { 'B5', 'C5', 'D5','E5', 'F5','G5'};
    cond3_3 = { 'B6', 'C6', 'D6','E6', 'F6','G6'};
    cond3_2 = { 'B7', 'C7', 'D7','E7', 'F7','G7'};
    cond2_3 = { 'B8', 'C8', 'D8','E8', 'F8','G8'};
    cond2_2 = { 'B9', 'C9', 'D9','E9', 'F9','G9'};
    condUT = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    cond2_1 = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
    if contains(traj(i).well, cond4_4 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 4;
        traj(i).prevdose = [75,75]; % note if more than two treatments this could be a vector
        traj(i).doseints = [4,4]; % same with this 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000; % this is from the number of cells intended to be seeded (from well label)
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond4_3 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 3;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [4,3];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond3_4 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 4;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [3,4];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond4_2 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 2;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [4,2];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond3_3 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 3;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [3,3];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond3_2 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 2;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [3,2];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond2_3 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 3;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [2,3];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond2_2 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 2;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [2,2];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end

    if contains(traj(i).well, cond2_1 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 1;
        traj(i).prevdose = [75,75]; 
        traj(i).doseints = [2, 1];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
        traj(i).seed = 2000;
    end
    
    if contains(traj(i).well, condUT )
        traj(i).dose = 75;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
end

%% Load in 3rd data set (1818: single dose)
for i = sz(2,3)+1:sz(3,3)-6 % size matrix, first row, second column
    k=i-sz(2,3);
    traj(i).time = N3(1:end,1);
    traj(i).rawN = N3(1:end,k+1);
    traj(i).date = '8-16-18'; % get this from original excel file
    traj(i).welllabel = T3(1, k+1);
    welllabel = string(T3(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 73;
end

for i = sz(2,3)+55:sz(3,3) % size matrix, first row, second column
    k=i-sz(2,3);
    traj(i).time = N3(1:end-1,1);
    traj(i).rawN = N3(1:end-1,k+1);
    traj(i).date = '8-16-18'; % get this from original excel file
    traj(i).welllabel = T3(1, k+1);
    welllabel = string(T3(1, k+1));
    wellfull = extractAfter(welllabel,"well "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = [];
    traj(i).tdose = 73;
end

  %%  Single Dose

    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
 for i = sz(2,3)+1:sz(3,3) % size matrix, first row, second column
    k=i-sz(2,3);
    cond300 = { 'B2', 'C2', 'D2','E2', 'F2','G2'};
    cond150 = { 'B3', 'C3', 'D3','E3', 'F3','G3'};
    cond125 = { 'B4', 'C4', 'D4','E4', 'F4','G4'};
    cond100 = { 'B5', 'C5', 'D5','E5', 'F5','G5'};
    cond75 = { 'B6', 'C6', 'D6','E6', 'F6','G6'};
    cond50 = { 'B7', 'C7', 'D7','E7', 'F7','G7'};
    cond35 = { 'B8', 'C8', 'D8','E8', 'F8','G8'};
    cond20 = { 'B9', 'C9', 'D9','E9', 'F9','G9'};
    cond10 = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    condUT = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
    if contains(traj(i).well, cond300 )
        traj(i).dose = 300;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; % note if more than two treatments this could be a vector
        traj(i).doseints = []; % same with this 
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000; % this is from the number of cells intended to be seeded (from well label)
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond150 )
        traj(i).dose = 150;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond125 )
        traj(i).dose = 125;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints =[];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond100 )
        traj(i).dose = 100;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond75 )
        traj(i).dose = 75;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond50 )
        traj(i).dose = 50;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond35 )
        traj(i).dose = 35;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
    if contains(traj(i).well, cond20 )
        traj(i).dose = 20;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end

    if contains(traj(i).well, cond10 )
        traj(i).dose = 10;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
        traj(i).flag = 0;
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 0;
        traj(i).dosenum = 0;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = [];
        traj(i).numdoses = 0;
        traj(i).seed = 2000;
=======
for i = 1:length(raw)
    sz(i,1:2) = size(raw(i).N);
    sz(i,2) = sz(i,2) - 1;
end

sz(:,3) = cumsum(sz(:,2));
%% Structure using fors and ifs
for i = 1:length(raw) %each data set runs
    if i == 1
        x = 1;
    else
        x = sz(i-1,3)+1;
    end
    for j = x:sz(i,3) % set bound in traj for data set
        % find well #, specific to the data set, and load data
        k = j + 1 - x;
        traj(j).time = raw(i).N(1:end,1);
        traj(j).rawN = raw(i).N(1:end,k+1);
        traj(j).N0true = traj(j).rawN(1);
        % load and mess with raw well labels to get actual well
        traj(j).welllabel =  raw(i).T(1,k+1);
        wellhalf = extractAfter(traj(j).welllabel,"(");
        wellact = extractAfter(wellhalf,"(");
        traj(j).well = strtok(wellact, {')'});
        traj(j).column = strtok(traj(j).well, {'B','C','D','E','F','G'});
        % subject to change; currently all data follows these
        traj(j).celltype = 'MCF-7';
        traj(j).drug = 'dox';
        traj(j).doseduration = 24;
        % individual data set parameters
        if i == 1       % GH1831: 2nd dose varying interval
            traj(j).dose = 75;
            traj(j).dosenum = 2;
            traj(j).numdoses = 2;
            traj(j).prevdose = 75;
            traj(j).tdose = 60;
            traj(j).date = '12-19-18';
            traj(j).seed = 2000;
            if ismember(traj(j).column, '2')
                traj(j).WPT = 8;
                traj(j).doseints = 8;
            elseif ismember(traj(j).column, '3')
                traj(j).WPT = 7;
                traj(j).doseints = 7;
            elseif ismember(traj(j).column, '4')
                traj(j).WPT = 6;
                traj(j).doseints = 6;
            elseif ismember(traj(j).column, '5')
                traj(j).WPT = 5;
                traj(j).doseints = 5;
            elseif ismember(traj(j).column, '6')
                traj(j).WPT = 4;
                traj(j).doseints = 4;
            elseif ismember(traj(j).column, '7')
                traj(j).WPT = 3;
                traj(j).doseints = 3;
            elseif ismember(traj(j).column, '8')
                traj(j).WPT = 2;
                traj(j).doseints = 1;
            elseif ismember(traj(j).column, '9')
                traj(j).WPT = 1;
                traj(j).doseints = 1;
            elseif ismember(traj(j).column, '10')
                traj(j).prevdose = [75, 100, 75, 75, 75]; 
                traj(j).doseints = [2, 4, 3, 4, 8];
                traj(j).dosenum = 6;
                traj(j).numdoses = 6;
                traj(j).WPT = 8
            elseif ismember(traj(j).column, '11')
                traj(j).dosenum = 1;
                traj(j).numdoses = 1;
                traj(j).WPT = [];
                traj(j).prevdose = [];
                traj(j).doseints = [];
            end
        elseif i == 2   %GH1830; 3rd dose varying interval
            traj(j).date = '12-19-18';
            traj(j).tdose = 69;
            traj(j).seed = 2000;
            traj(j).dosenum = 3;
            traj(j).numdoses = 3;
            traj(j).prevdose = [75,75];
            traj(j).dose = 75;
            if ismember(traj(j).column, '2')
                traj(j).doseints = [4,4]
                traj(j).WPT = 4;
            elseif ismember(traj(j).column, '3')
                traj(j).doseints = [4,3]
                traj(j).WPT = 3;
            elseif ismember(traj(j).column, '4')
                traj(j).doseints = [3,4]
                traj(j).WPT = 4;
            elseif ismember(traj(j).column, '5')
                traj(j).doseints = [4,2]
                traj(j).WPT = 2;
            elseif ismember(traj(j).column, '6')
                traj(j).doseints = [3,3]
                traj(j).WPT = 3;
            elseif ismember(traj(j).column, '7')
                traj(j).doseints = [3,2]
                traj(j).WPT = 2;
            elseif ismember(traj(j).column, '8')
                traj(j).doseints = [2,3]
                traj(j).WPT = 3;
            elseif ismember(traj(j).column, '9')
                traj(j).doseints = [2,2]
                traj(j).WPT = 2;
            elseif ismember(traj(j).column, '10')
                traj(j).doseints = [2,1]
                traj(j).WPT = 4;
            elseif ismember(traj(j).column, '11')
                traj(j).dosenum = 1;
                traj(j).numdoses = 1;
                traj(j).WPT = [];
                traj(j).prevdose = [];
                traj(j).doseints = [];
            end
        elseif i == 3 %GH1818 single dose
            traj(j).date = '8=16=18';
            traj(j).tdose = 73;
            traj(j).seed = 2000;
            traj(j).dosenum = 1;
            traj(j).prevdose = [];
            traj(j).WPT = [];
            traj(j).numdoses = 1;
            traj(j).doseints = [];
            if ismember(traj(j).column, '2')
                traj(j).dose = 300;
            elseif ismember(traj(j).column, '3')
                traj(j).dose = 150;
            elseif ismember(traj(j).column, '4')
                traj(j).dose = 125;
            elseif ismember(traj(j).column, '5')
                traj(j).dose = 100;
            elseif ismember(traj(j).column, '6')
                traj(j).dose = 75;
            elseif ismember(traj(j).column, '7')
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '8')
                traj(j).dose = 35;
            elseif ismember(traj(j).column, '9')
                traj(j).dose = 20;
            elseif ismember(traj(j).column, '10')
                traj(j).dose = 10;
            elseif ismember(traj(j).column, '11')   
                traj(j).dose = [];
                traj(j).dosenum = [];
                traj(j).numdoses = 0;
            end
        elseif i == 4 %GH1832 3rd dose varying dosage
            traj(j).date = '12-18-18'
            traj(j).tdose = 69;
            traj(j).seed = 2000;
            traj(j).dosenum = 3;
            traj(j).numdoses = 3;
            traj(j).WPT = 4;
            traj(j).prevdose = [75,75];
            traj(j).doseints = [4,4];
            if ismember(traj(j).column, '2')
                traj(j).dose = 300;
            elseif ismember(traj(j).column, '3')
                traj(j).dose = 150;
            elseif ismember(traj(j).column, '4')
                traj(j).dose = 125;
            elseif ismember(traj(j).column, '5')
                traj(j).dose = 100;
            elseif ismember(traj(j).column, '6')
                traj(j).dose = 75;
            elseif ismember(traj(j).column, '7')
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '8')
                traj(j).dose = 35;
            elseif ismember(traj(j).column, '9')
                traj(j).dose = 20;
            elseif ismember(traj(j).column, '10')
                traj(j).dose = 10;
            elseif ismember(traj(j).column, '11')
                traj(j).dose = 0;
            end
        elseif i == 5 %GH1820: 2nd/3rd dose varying interval
            traj(j).date = '9-15-18';
            traj(j).tdose = 63;
            traj(j).seed = 3000;
            traj(j).dose = 75;
            traj(j).dosenum = 2;
            traj(j).numdoses = 2;
            traj(j).prevdose = 75;
            if ismember(traj(j).column, '2')
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).WPT = 4;
                traj(j).prevdose = [75,100]
                traj(j).doseints = [6,4];
            elseif ismember(traj(j).column, '3')
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).WPT = 1;
                traj(j).prevdose = [75,50];
                traj(j).doseints = [5,1];
            elseif ismember(traj(j).column, '4')
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).WPT = 1;
                traj(j).prevdose = [75,35];
                traj(j).doseints = [4,1];
            elseif ismember(traj(j).column, '5')
                traj(j).WPT = 6;
                traj(j).doseints = 6;
            elseif ismember(traj(j).column, '6')
                traj(j).WPT = 5;
                traj(j).doseints = 5;
            elseif ismember(traj(j).column, '7')
                traj(j).WPT = 4;
                traj(j).doseints = 4;
            elseif ismember(traj(j).column, '8')
                traj(j).WPT = 3;
                traj(j).doseints = 3;
            elseif ismember(traj(j).column, '9')
                traj(j).WPT = 2;
                traj(j).doseints = 2;
            elseif ismember(traj(j).column, '10')
                traj(j).WPT = 1;
                traj(j).doseints = 1;
            elseif ismember(traj(j).column, '11') 
                traj(j).dosenum = 1;
                traj(j).numdoses = 1;
                traj(j).WPT = [];
                traj(j).prevdose = [];
                traj(j).doseints = [];
            end    
        elseif i == 6 %GH1829: 2nd dose varying dosage
            traj(j).date = '12-11-18';
            traj(j).tdose = 68.38;
            traj(j).seed = 2000;
            traj(j).dosenum = 2;
            traj(j).numdoses = 2;
            traj(j).WPT = 4;
            traj(j).prevdose = 75;
            traj(j).doseints = 4;
            if ismember(traj(j).column, '2')
                traj(j).dose = 300;
            elseif ismember(traj(j).column, '3')
                traj(j).dose = 150;
            elseif ismember(traj(j).column, '4')
                traj(j).dose = 125;
            elseif ismember(traj(j).column, '5')
                traj(j).dose = 100;
            elseif ismember(traj(j).column, '6')
                traj(j).dose = 75;
            elseif ismember(traj(j).column, '7')
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '8')
                traj(j).dose = 35;
            elseif ismember(traj(j).column, '9')
                traj(j).dose = 20;
            elseif ismember(traj(j).column, '10')
                traj(j).dose = 10;
            elseif ismember(traj(j).column, '11')   
                traj(j).dose = 0;
            end    
        elseif i == 7 %GH1839 3rd dose varying dosage
            traj(j).date = '8-23-18';
            traj(j).tdose = 72;
            traj(j).seed = 2000;
            traj(j).dosenum = 3;
            traj(j).numdoses = 3;
            traj(j).WPT = 2;
            traj(j).prevdose = [75,100];
            traj(j).doseints = [2,2];
            if ismember(traj(j).column, '2')
                traj(j).dose = 300;
            elseif ismember(traj(j).column, '3')
                traj(j).dose = 150;
            elseif ismember(traj(j).column, '4')
                traj(j).dose = 125;
            elseif ismember(traj(j).column, '5')
                traj(j).dose = 100;
            elseif ismember(traj(j).column, '6')
                traj(j).dose = 75;
            elseif ismember(traj(j).column, '7')
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '8')
                traj(j).dose = 35;
            elseif ismember(traj(j).column, '9')
                traj(j).dose = 20;
            elseif ismember(traj(j).column, '10')
                traj(j).dose = 10;
            elseif ismember(traj(j).column, '11')   
                traj(j).dose = 0;
            end    
        elseif i == 8 %GH1817: 2nd dose varying dosage
            traj(j).date = '8-19-18';
            traj(j).tdose = 68;
            traj(j).seed = 2000;
            traj(j).dosenum = 2;
            traj(j).numdoses = 2;
            traj(j).WPT = 2;
            traj(j).prevdose = 75;
            traj(j).doseints = 2;
            if ismember(traj(j).column, '2')
                traj(j).dose = 300;
            elseif ismember(traj(j).column, '3')
                traj(j).dose = 150;
            elseif ismember(traj(j).column, '4')
                traj(j).dose = 125;
            elseif ismember(traj(j).column, '5')
                traj(j).dose = 100;
            elseif ismember(traj(j).column, '6')
                traj(j).dose = 75;
            elseif ismember(traj(j).column, '7')
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '8')
                traj(j).dose = 35;
            elseif ismember(traj(j).column, '9')
                traj(j).dose = 20;
            elseif ismember(traj(j).column, '10')
                traj(j).dose = 10;
            elseif ismember(traj(j).column, '11')   
                traj(j).dose = 0;
            end
        elseif i == 9 %GH1815: 2nd dose varying initial dosage
            traj(j).date = '7-12-18';
            traj(j).tdose = 76.55
            traj(j).seed = 2000;
            traj(j).numdoses = 2;
            traj(j).dosenum = 2;
            traj(j).doseints = 2;
            traj(j).WPT = 2;
            traj(j).dose = 100;
            if ismember(traj(j).column, '2')
                traj(j).dosenum = 1;
                traj(j).numdoses = 1;
                traj(j).WPT = [];
                traj(j).doseints = [];
                traj(j).prevdose = [];
            elseif ismember(traj(j).column, '3')
                traj(j).prevdose = 10;
            elseif ismember(traj(j).column, '4')
                traj(j).prevdose = 20;
            elseif ismember(traj(j).column, '5')
                traj(j).prevdose = 35;
            elseif ismember(traj(j).column, '6')
                traj(j).prevdose = 50;
            elseif ismember(traj(j).column, '7')
                traj(j).prevdose = 75;
            elseif ismember(traj(j).column, '8')
                traj(j).prevdose = 100;
            elseif ismember(traj(j).column, '9')
                traj(j).prevdose = 125;
            elseif ismember(traj(j).column, '10')
                traj(j).prevdose = 150;
            elseif ismember(traj(j).column, '11')    
                traj(j).prevdose = 300;
            end
        elseif i == 10 %GH1835: 2nd dose varying interval
            traj(j).date = '12-18-18';
            traj(j).tdose = 60;
            traj(j).celltype = 'BT474';
            traj(j).seed = 2000;
            traj(j).dosenum = 2;
            traj(j).numdoses = 2;
            traj(j).dose = 35;
            traj(j).prevdose = 35;
            if ismember(traj(j).column, '2')
                traj(j).WPT = 7;
                traj(j).doseints = 7;
            elseif ismember(traj(j).column, '3')
                traj(j).WPT = 6;
                traj(j).doseints = 6;
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '4')
                traj(j).WPT = 6;
                traj(j).doseints = 6;
            elseif ismember(traj(j).column, '5') 
                traj(j).WPT = 5;
                traj(j).doseints = 5;
                traj(j).dose = 50;
            elseif ismember(traj(j).column, '6')
                traj(j).WPT = 5;
                traj(j).doseints = 5;
            elseif ismember(traj(j).column, '7')
                traj(j).WPT = 4;
                traj(j).doseints = 4;
            elseif ismember(traj(j).column, '8')
                traj(j).WPT = 3;
                traj(j).doseints = 3;
            elseif ismember(traj(j).column, '9')
                traj(j).WPT = 2;
                traj(j).doseints = 2;
            elseif ismember(traj(j).column, '10')
                traj(j).WPT = 1;
                traj(j).doseints = 1;
            elseif ismember(traj(j).column, '11')
                traj(j).WPT = []
                traj(j).doseints = [];
                traj(j).dosenum = 1;
                traj(j).numdoses = 1;
                traj(j).prevdose = [];
            end
        end
        traj(j).accdose = traj(j).dose + sum(traj(j).prevdose);
>>>>>>> d87b9cf12f9dc23e92ac5061858377a2ac9ee47b
    end
end
    


%% Add color by WPT 
for j = 1:length(traj)
    if ~isempty(traj(j).WPT)
    WPT(j) = traj(j).WPT;
    end
end
colorsets = varycolor(length(unique(WPT))+2);
uniqWPT= unique(WPT);

for i = 1:length(traj)
    traj(i).color = [];
    for j = 1:length(uniqWPT)
        if traj(i).WPT==uniqWPT(j)
            %used to be if flag = 0, none are flagged now
            traj(i).color =colorsets(j,:);
        end
    end
    if isempty(traj(i).WPT)
        traj(i).color = colorsets(end,:); % make untreated control black
    end
end

%% Add color by 1st dose
for j = 1:length(traj)
    if traj(j).numdoses==1
    if ~isempty(traj(j).dose)
    dose(j) = traj(j).dose;
    end
    end
end
colorsets = varycolor(length(unique(dose))+1);
uniqdose= unique(dose);

for i = 1:length(traj)
   
    for j = 1:length(uniqdose)
        if traj(i).numdoses ==1
        if traj(i).dose==uniqdose(j)
            %used to be if flag = 0, now none are flagged
            traj(i).color =colorsets(j,:);
        end
        end
    end
    if isempty(traj(i).dose)
        traj(i).color = colorsets(end,:); % make untreated control black
    end
end
%% Save the raw data structure data sets
% This saves the traj structure just containing raw data as trajraw.mat
save('../out/trajraw.mat', 'traj')