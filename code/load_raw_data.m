% This script loads in data sets from Grant's various experimental
% conditions of drug response. The goal of this script is just to capture
% the input variables and data for each well in the data set.

close all; clear all; clc

% load in adjusted data set (remove top part and record relevant info
% below)
[N1, T1] =xlsread('../data/GH1831_MCF7_PredModel_2nd_Dose_Varying_Interval_v2.xls');
[N2, T2] = xlsread('../data/GH1830_MCF7_PredModel_3rd_Dose_Vary_Interval_v2.xls');
[N3, T3] = xlsread('../data/GH1818_MCF7_PredModel_1st_Dose_Analysis_v6.xls');
[N4, T4] =xlsread('../data/GH1832_MCF7_PredModel_3rd_Dose_4-4_Intervals_v1.xls');
[N5, T5] =xlsread('../data/GH1820_MCF7_PredModel_2nd_Dose_Varying_Interval_v3.xls');
[N6, T6] =xlsread('../data/GH1829_MCF7_PredModel_2nd_Dose_4_Week_Interval_v1.xls');
[N7, T7] =xlsread('../data/GH1819_MCF7_PredModel_3rd_Treatment_Analysis_v2.xls');
[N8, T8] =xlsread('../data/GH1817_MCF7_PredModel_2nd_Dose_Analysis_v3.xls');
[N9, T9] =xlsread('../data/GH1815_MCF7_PredModel_Effect_of_1st_Dose_Analysis_v2.xls');

% this reads in the numeric values in an excel sheet (N1), and the text
% values (T1)
% Make a "size" matrix that keeps track of the number of wells and time
% points of each data set loaded. This will be used to add to the structure
% once more data sets are loaded in
% 1st column in number of time points, 2nd column number of wells in that
% data set, third column is cumulative number of wells
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
<<<<<<< HEAD
  %%  Single Dose
=======
>>>>>>> e37c770e18c9ce972504cb71c27614e9a06373aa
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
    end
end
%% Load in fourth data set (1832: 3rd dose with different intervals)

for i = sz(3,3)+1:sz(4,3) % size matrix, first row, second column
    k=i-sz(3,3);
    traj(i).time = N4(1:end,1);
    traj(i).rawN = N4(1:end,k+1);
    traj(i).date = '12-18-18'; % get this from original excel file
    traj(i).welllabel = T4(1, k+1);
    welllabel = string(T4(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the
    well = strtok(wellfull, {'(', ')'});
    if k > 54
        wellfull = extractAfter(welllabel,"well "); % here want to extract the
        well = strtok(wellfull, {'(', ')'});
        traj(i).well = well;
    end
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
    
    % All wells in this plate have identical initial conditions
    traj(i).dosenum = 3;
    traj(i).WPT = 4;
    traj(i).prevdose = [75,75]; % note if more than two treatments this could be a vector
    traj(i).doseints = [4,4]; % same with this 
    traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    traj(i).numdoses = 3;
    traj(i).seed = 2000; % this is from the number of cells intended to be seeded (from well label)
    
    if contains(traj(i).well, cond300 )
        traj(i).dose = 300;
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond150 )
        traj(i).dose = 150;
    end
    if contains(traj(i).well, cond125 )
        traj(i).dose = 125;
    end
    if contains(traj(i).well, cond100 )
        traj(i).dose = 100;
    end
    if contains(traj(i).well, cond75 )
        traj(i).dose = 75;
    end
    if contains(traj(i).well, cond50 )
        traj(i).dose = 50;
    end
    if contains(traj(i).well, cond35 )
        traj(i).dose = 35;
    end
    if contains(traj(i).well, cond20 )
        traj(i).dose = 20;
    end
    if contains(traj(i).well, cond10 )
        traj(i).dose = 10;
    end  
    if contains(traj(i).well, condUT )
        traj(i).dose = 0;
        traj(i).dosenum = []; % override "initial" settings that applied to the other wells
        traj(i).numdoses = 2;
        
    end
end

%% Fifth data set (1820: second dose varying interval)
for i = sz(4,3) + 1:sz(5,3) % size matrix, first row, second column
    k = i - sz(4,3);
    traj(i).time = N5(1:end,1);
    traj(i).rawN = N5(1:end,k+1);
    traj(i).date = '9-5-18'; % get this from original excel file
    traj(i).welllabel = T5(1, k+1);
    welllabel = string(T5(1, k+1));
    wellfull = extractAfter(welllabel,"well "); % here want to extract the
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 63;
    traj(i).seed = 3000; % this is from the number of cells intended to be seeded (from well label)
 
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
    cond64 = { 'B2', 'C2', 'D2','E2', 'F2','G2'};
    cond51 = { 'B3', 'C3', 'D3','E3', 'F3','G3'};
    cond41 = { 'B4', 'C4', 'D4','E4', 'F4','G4'};
    cond6 = { 'B5', 'C5', 'D5','E5', 'F5','G5'};
    cond5 = { 'B6', 'C6', 'D6','E6', 'F6','G6'};
    cond4 = { 'B7', 'C7', 'D7','E7', 'F7','G7'};
    cond3 = { 'B8', 'C8', 'D8','E8', 'F8','G8'};
    cond2 = { 'B9', 'C9', 'D9','E9', 'F9','G9'};
    cond1 = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    condUT = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
    if contains(traj(i).well, cond64 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 4;
        traj(i).prevdose = [75,100]; % note if more than two treatments this could be a vector
        traj(i).doseints = [6,4]; % same with this 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond51)
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 1;
        traj(i).prevdose = [75,50]; 
        traj(i).doseints = [5,1];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
    end
    if contains(traj(i).well, cond41 )
        traj(i).dose = 75;
        traj(i).dosenum = 3;
        traj(i).WPT = 1;
        traj(i).prevdose = [75,35]; 
        traj(i).doseints = [4,1];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 3;
    end
    if contains(traj(i).well, cond6 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 6;
        traj(i).prevdose = 75; 
        traj(i).doseints = 6;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, cond5 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 5;
        traj(i).prevdose = 75; 
        traj(i).doseints = 5;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, cond4 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 4;
        traj(i).prevdose = 75; 
        traj(i).doseints = 4;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, cond3 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 3;
        traj(i).prevdose = 75; 
        traj(i).doseints = 3;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, cond2 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 2;
        traj(i).prevdose = 75; 
        traj(i).doseints = 2;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, cond1 )
        traj(i).dose = 75;
        traj(i).dosenum = 2;
        traj(i).WPT = 1;
        traj(i).prevdose = 75; 
        traj(i).doseints = 1;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 75;
        traj(i).dosenum = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; 
        traj(i).doseints = [];
        traj(i).accdose = traj(i).dose;
        traj(i).numdoses = 1;
    end
end
%% Load in 6th data set (1829: second dose varying dosage)
    % initial seeding of 2k is only assumed
for i = sz(5,3)+1:sz(6,3) - 6% size matrix, first row, second column
    k=i-sz(5,3);
    traj(i).time = N6(1:end,1);
    traj(i).rawN = N6(1:end,k+1);
    traj(i).date = '12-11-18'; % get this from original excel file
    traj(i).welllabel = T6(1, k+1);
    welllabel = string(T6(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 68.38;
    traj(i).WPT = 4;
    traj(i).dosenum = 2;
    traj(i).prevdose = 75;
    traj(i).doseints = 4;
    traj(i).numdoses = 2;
    traj(i).seed = 2000;
end
for i = sz(5,3)+55:sz(6,3) % size matrix, first row, second column
    k=i-sz(5,3);
    traj(i).time = N6(1:end-1,1);
    traj(i).rawN = N6(1:end-1,k+1);
    traj(i).date = '12-11-18'; % get this from original excel file
    traj(i).welllabel = T6(1, k+1);
    welllabel = string(T6(1, k+1));
    wellfull = extractAfter(welllabel,"well "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = [];
    traj(i).tdose = [];
end
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
for i = sz(5,3)+1:sz(6,3)
    k = i-sz(5,3);
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
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond150 )
        traj(i).dose = 150;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond125 )
        traj(i).dose = 125;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond100 )
        traj(i).dose = 100;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond75 )
        traj(i).dose = 75;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond50 )
        traj(i).dose = 50;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond35 )
        traj(i).dose = 35;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond20 )
        traj(i).dose = 20;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond10 )
        traj(i).dose = 10;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 0;
        traj(i).dosenum = [];
        traj(i).WPT = 4;
        traj(i).prevdose = 75; 
        traj(i).doseints = 4;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
end
%% Load in 7th data set (1819: third dose varying dosage)
    % initial seeding of 2k is only assumed
for i = sz(6,3)+1:sz(7,3) - 6% size matrix, first row, second column
    k=i-sz(6,3);
    traj(i).time = N7(1:end,1);
    traj(i).rawN = N7(1:end,k+1);
    traj(i).date = '8-23-18'; % get this from original excel file
    traj(i).welllabel = T7(1, k+1);
    welllabel = string(T7(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 72;
    traj(i).WPT = 2;
    traj(i).dosenum = 3;
    traj(i).prevdose =[75,100];
    traj(i).doseints = [2,2];
    traj(i).numdoses = 3;
    traj(i).seed = 2000;
end
for i = sz(6,3)+55:sz(7,3) % size matrix, first row, second column
    k=i-sz(6,3);
    traj(i).time = N7(1:end-1,1);
    traj(i).rawN = N7(1:end-1,k+1);
    traj(i).date = '8-23-18'; % get this from original excel file
    traj(i).welllabel = T7(1, k+1);
    welllabel = string(T7(1, k+1));
    wellfull = extractAfter(welllabel,"well "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = [];
    traj(i).tdose = [];
end
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
for i = sz(6,3)+1:sz(7,3)
    k = i-sz(6,3);
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
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond150 )
        traj(i).dose = 150;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond125 )
        traj(i).dose = 125;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond100 )
        traj(i).dose = 100;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond75 )
        traj(i).dose = 75;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond50 )
        traj(i).dose = 50;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond35 )
        traj(i).dose = 35;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond20 )
        traj(i).dose = 20;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond10 )
        traj(i).dose = 10;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 0;
        traj(i).dosenum = [];
        traj(i).WPT = 2;
        traj(i).prevdose = [75,100]; 
        traj(i).doseints = [2,2];
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 2;
        traj(i).seed = 2000;
    end
end
%% Load in 8th data set (1817: second dose varying dosage, 2 WPT)
    % initial seeding of 2k is only assumed
for i = sz(7,3)+1:sz(8,3) - 6% size matrix, first row, second column
    k=i-sz(7,3);
    traj(i).time = N8(1:end,1);
    traj(i).rawN = N8(1:end,k+1);
    traj(i).date = '8-9-18'; % get this from original excel file
    traj(i).welllabel = T8(1, k+1);
    welllabel = string(T8(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 68;
    traj(i).WPT = 2;
    traj(i).dosenum = 2;
    traj(i).prevdose = 75;
    traj(i).doseints = 2;
    traj(i).numdoses = 2;
    traj(i).seed = 2000;
end
for i = sz(7,3)+55:sz(8,3) % size matrix, first row, second column
    k=i-sz(7,3);
    traj(i).time = N8(1:end-1,1);
    traj(i).rawN = N8(1:end-1,k+1);
    traj(i).date = '8-9-18'; % get this from original excel file
    traj(i).welllabel = T8(1, k+1);
    welllabel = string(T8(1, k+1));
    wellfull = extractAfter(welllabel,"well "); % here want to extract the well
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = [];
    traj(i).tdose = [];
end
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
for i = sz(7,3)+1:sz(8,3)
    k = i-sz(7,3);
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
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond150 )
        traj(i).dose = 150;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond125 )
        traj(i).dose = 125;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond100 )
        traj(i).dose = 100;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond75 )
        traj(i).dose = 75;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond50 )
        traj(i).dose = 50;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond35 )
        traj(i).dose = 35;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond20 )
        traj(i).dose = 20;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond10 )
        traj(i).dose = 10;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, condUT )
        traj(i).dose = 0;
        traj(i).dosenum = [];
        traj(i).WPT = 2;
        traj(i).prevdose = 75; 
        traj(i).doseints = 2;
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
        traj(i).numdoses = 1;
        traj(i).seed = 2000;
    end
end
%% 9th data set (1815: 2nd dose with different initial dose)

for i = sz(8,3)+1:sz(9,3) % size matrix, first row, second column
    k=i-sz(8,3);
    traj(i).time = N9(1:end,1);
    traj(i).rawN = N9(1:end,k+1);
    traj(i).date = '7-12-18'; % get this from original excel file
    traj(i).welllabel = T9(1, k+1);
    welllabel = string(T9(1, k+1));
    wellfull = extractAfter(welllabel,"nM "); % here want to extract the
    if k > 6
        welllabel = wellfull;
        wellfull = extractAfter(welllabel,"nM "); % there are two "nM"s
    end
    well = strtok(wellfull, {'(', ')'});
    traj(i).well = well;
    traj(i).flag =0;
    traj(i).celltype = 'MCF-7';
    traj(i).drug = 'dox';
    traj(i).N0true= traj(i).rawN(1);
    traj(i).doseduration = 24;
    traj(i).tdose = 76.55;
    traj(i).dose = 100;
    traj(i).dosenum = 2;
    traj(i).numdoses = 2;
    traj(i).WPT = 2;
    traj(i).doseints = 2;
    traj(i).seed = 2000;
 
    
    % Need to write in a little "key" for this data set to match wells to
    % their experimental condition. This is the relatively tedious part....
    % Need to go through original excel sheet and code in which wells
    % correspond to which conditions
    condUT = { 'B2', 'C2', 'D2','E2', 'F2','G2'};
    cond10 = { 'B3', 'C3', 'D3','E3', 'F3','G3'};
    cond20 = { 'B4', 'C4', 'D4','E4', 'F4','G4'};
    cond35 = { 'B5', 'C5', 'D5','E5', 'F5','G5'};
    cond50 = { 'B6', 'C6', 'D6','E6', 'F6','G6'};
    cond75 = { 'B7', 'C7', 'D7','E7', 'F7','G7'};
    cond100 = { 'B8', 'C8', 'D8','E8', 'F8','G8'};
    cond125 = { 'B9', 'C9', 'D9','E9', 'F9','G9'};
    cond150 = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    cond300 = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
    if contains(traj(i).well, condUT )
        traj(i).dosenum = 1;
        traj(i).numdoses = 1;
        traj(i).WPT = [];
        traj(i).prevdose = []; % note if more than two treatments this could be a vector
        traj(i).doseints = []; % same with this 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    % Repeat this for all other conditions (usually columns of a plate)
    if contains(traj(i).well, cond10 )
        traj(i).prevdose = 10; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond20 )
        traj(i).prevdose = 20; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond35)
        traj(i).prevdose = 35; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond50)
        traj(i).prevdose = 50; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond75)
        traj(i).prevdose = 75; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond100)
        traj(i).prevdose = 100; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond125 )
        traj(i).prevdose = 125; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond150)
        traj(i).prevdose = 150; 
        traj(i).accdose = traj(i).dose + sum(traj(i).prevdose);
    end
    if contains(traj(i).well, cond300 )
        traj(i).prevdose = 300; 
        traj(i).accdose = traj(i).dose;
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
            if traj(i).flag ==0
            traj(i).color =colorsets(j,:);
            end
        end
    end
    if isempty(traj(i).WPT)
        traj(i).color = colorsets(end,:); % make untreated control black
    end
    if traj(i).flag==1
        traj(i).color = colorsets(end-1, :);
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
            if traj(i).flag ==0
            traj(i).color =colorsets(j,:);
            end
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