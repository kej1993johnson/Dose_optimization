%function [profiles] = profile_likelihood(params, tsamp, N0, N,mudatavec, vardatavec, modelcode, factor, numpoints)
%fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
%% Play with this to set different parameters

% What if we set rr = xx from Daylins AA 161 lineage?
psetID = [1, 3, 4, 6]% phi0, carcapN, carcapphi, rr
pfitID = [2, 5, 7, 8]; % corresponds to rs, alpha, ds, dr
% Forces rr to be much higher than it is fitting as
rr_to_pop_ratio= 0.0317/0.0392;
rrdata = rr_to_pop_ratio*gtot
pset = [phi0, carcapNf, carcapphi, rrdata];
%%
nfitparams = length(pfitID);
        profile = [];
factor = 0.2;
numpoints = 10;
params = pbest; % rs, alpha, zr, ds, zd];


%%
for k = 1:nfitparams
profindex = pfitID(k); % PROFILE the kth fit parameter 
    
    profrangeDown = linspace((params(k)), (params(k)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(k)), (params(k)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange);
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    profile = [];
    
    for m = 1:length(profrange)
        [m] %track progress
        currp = profrange(m);
        % Change pfitID and psetID so that pfit is all params except
        % currently profiled param and psetID includes currently profiled
        % param
        [psetIDcurr, ordp] = sort([psetID, profindex]);
        ikeep = pfitID~=profindex; % identify all indices corresponding to
        % parameters not being profiled
        pfitIDcurr = pfitID(ikeep);
        thetacurr = pbest(ikeep);
        pboundscurr = pbounds(ikeep, :);
        pcomb = [pset, currp];
        psetcurr = pcomb(ordp);
        [paramstemp,~, ~, fvaltemp, ~, ~]=fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitIDcurr, psetIDcurr, thetacurr, psetcurr, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pboundscurr);
        
        %[fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly
      
        currfvals = [currfvals; fvaltemp];
        currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
    
    profile = horzcat(profrange, real(currfvals));
 


    profiles(:,:,k) = profile;

% 1st column is the parameter values that are "profiled"
% 2nd column is the negLL corresponding to each "profiled" parameter
% each dimension is for a parameter




 
end

%%
figure;
plot(profile(:,1), profile(:,2),'b-')
hold on
plot(pbest(2), negLL, 'r*')
plot([profile(1,1) profile(end,1)],[threshold threshold],'r--')
xlabel('profiled param')
ylabel('negLL')
set(gca,'FontSize',20,'LineWidth',1.5)
title(rs)
%% Plot the profiles
%plist = {'r_{s}', '\alpha', 'r_{ratio}', 'd_{s}', 'd_{ratio}'};
plist = {'\phi_{0}','r_{s}', 'r_{r}/r_{s} ratio', '\alpha', 'd_{s}', 'd_{r}/d_{s} ratio'};
threshold = chi2inv(0.95,length(pfitID))/2 + negLL;
figure;
for i = 1:length(pbest)
subplot(1,length(pbest),i)
plot(profiles(:,1, i), profiles(:,2, i),'b-', 'LineWidth', 1.8)
hold on
plot(pbest(i), negLL, 'r*')
plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
xlabel('profiled param')
ylabel('negLL')
set(gca,'FontSize',20,'LineWidth',1.5)
title(plist(i))
xlim([profiles(1,1,i), profiles(end,1, i)])
end
