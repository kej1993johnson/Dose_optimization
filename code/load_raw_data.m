% This script loads in data sets from Grant's various experimental
% conditions of drug response. The goal of this script is just to capture
% the input variables and data for each well in the data set.

close all; clear all; clc

% load in adjusted data set (remove top part and record relevant info
% below)
[N1, T1] =xlsread('../data/GH_MCF7_2nd_Dose_Vary_Interval_v1.xls');
[N2, T2] = xlsread('../data/GH1830_MCF7_PredModel_3rd_Dose_Vary_Interval_v2.xls');


% this reads in the numeric values in an excel sheet (N1), and the text
% values (T1)
% Make a "size" matrix that keeps track of the number of wells and time
% points of each data set loaded. This will be used to add to the structure
% once more data sets are loaded in
% 1st column in number of time points, 2nd column number of wells in that
% data set, third column is cumulative number of wells
sz(1,1:2) = size(N1);
sz(1,2)=sz(1,2)-1;
sz(:,3)= cumsum(sz(:,2));

sz(2,1:2) = size(N2);
sz(2,2) = sz(2,2)-1;
sz(:,3)= cumsum(sz(:,2));
%% Make structure
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
    traj(i).time = N1(1:end-1,1);
    traj(i).rawN = N1(1:end-1,i+1);
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
    traj(i).tdose = 60;
 
    
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
        traj(i).numdoses = 1;
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
%% Load in second data set ( 3rd dose with different intervals)

for i = sz(1,3)+1:sz(2,3) % size matrix, first row, second column
    k=i-sz(1,3);
    traj(i).time = N2(1:end-1,1);
    traj(i).rawN = N2(1:end-1,k+1);
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
    cond2_1 = { 'B10', 'C10', 'D10','E10', 'F10','G10'};
    condUT = { 'B11', 'C11', 'D11','E11', 'F11','G11'};
    
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
%% Save the raw data structure data sets
% This saves the traj structure just containing raw data as trajraw.mat
save('../out/trajraw.mat', 'traj')