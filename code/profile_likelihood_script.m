%function [profiles] = profile_likelihood(params, tsamp, N0, N,mudatavec, vardatavec, modelcode, factor, numpoints)
%fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
%% Play with this to set different parameters

% What if we set rr = xx from Daylins AA 161 lineage?
psetID = [1, 3, 4]% phi0, carcapN, carcapphi, rr
pfitID = [2, 5, 6, 7, 8]; % corresponds to rs, alpha, ds, dr
% Forces rr to be much higher than it is fitting as
rr_to_pop_ratio= 0.0317/0.0392;
rrdata = rr_to_pop_ratio*gtot
pset = [phi0, carcapNf, carcapphi];
%%
nfitparams = length(pfitID);
        profile = [];
factor0 = 0.6;
numpoints = 20;
params = pbest; % rs, alpha, zr, ds, zd];
threshold = chi2inv(0.95,length(pfitID))/2 + negLL;
profiles = [];
%% For 231s without fitting phi0
pboundsprof = [0, Inf; 0, Inf; 0, Inf; 0, Inf; 0, Inf];
%% For CLL
Ub = Uphi;
lengthvec = lengthvecN;
%% 231s with fitting phi0
pboundsprof = [0, Inf; 0, Inf; 0, Inf; 0, Inf; 0, Inf; 0, Inf];
%%

ikeep = [];
for k = 1:nfitparams
profindex = pfitID(k); % PROFILE the kth fit parameter 
    factor = factor0;
    if k == 5 % increase factor for ds and dr/ds ratio 
        factor = 2*factor0;
    end
    if k==6
        factor = 3*factor0;
    end
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
        pboundscurr = pboundsprof(ikeep, :);
        pcomb = [pset, currp];
        psetcurr = pcomb(ordp);
        [paramstemp,~, ~, fvaltemp, ~, ~]=fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitIDcurr, psetIDcurr, thetacurr, psetcurr, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pboundscurr);
        
        %[fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly
      
        currfvals = [currfvals; fvaltemp];
        %currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
        currparams = [currparams; [paramstemp(1:k-1),currp,paramstemp(k:end)]];
    end
    
    profile = horzcat(profrange, real(currfvals), currparams);
 


    profiles(:,:,k) = profile;
% each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS, negLL) value, and
        % then columns for each of the other parameter estimates.
    

ilo = [];
 ilo=find(profiles(:,2,k)<threshold);
     
        if ilo(end) ==numpoints*2 ||ilo(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(k,:) = [NaN, NaN];
        else
             CI(k,:) = [profiles(ilo(1),1,k), profiles(ilo(end),1,k)];
        end

end



%% Plot the profiles
plist = {'r_{s}', '\alpha', 'rr/rs', 'd_{s}', 'dr/ds'};
plist = {'\phi_{0}','r_{s}', 'r_{r}/r_{s} ratio', '\alpha', 'd_{s}', 'd_{r}/d_{s} ratio'};
threshold = chi2inv(0.95,length(pfitID))/2 + negLL;
figure;
for i = 1:length(pbest)
subplot(2,3,i)
plot(profiles(:,1, i), profiles(:,2, i),'b-', 'LineWidth', 1.8)
hold on
plot(pbest(i), negLL, 'r*')
plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
xlabel(plist(i))
ylabel('cost')
set(gca,'FontSize',16,'LineWidth',1.5)
%title(plist(i))
ylim([75 90])
xlim([profiles(1,1,i),profiles(end,1,i)])
%xlim([CI(i,1)-0.1*pbest(i),CI(i,2)+0.1*pbest(i)])
end
%% If your parameters are identifiable, save the CI so that we can use these later to compare data
% uncertainty to model uncertainty
save('../out/CIpbest.mat', 'CI')
%% Plot parameter relationships
for i = 1:length(pbest)
    figure;
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    plot(profiles(:,1,i),profiles(:,3:end,i),'LineWidth',2)
    plot(pbest(i),pbest,'r*')
    xlabel(plist{i})
    ylabel('Estimated Parameter Value')
    legend(plist)
    legend boxoff
    title(plist(i))
end
%% Is there another good way to plot the parameters and their CI?
