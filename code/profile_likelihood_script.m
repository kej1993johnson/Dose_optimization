%function [profiles] = profile_likelihood(params, tsamp, N0, N,mudatavec, vardatavec, modelcode, factor, numpoints)
%fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);

nfitparams = length(pbest);
        profile = [];
factor = 0.1;
numpoints = 10;
params = pbest; % rs, alpha, zr, ds, zd];

%%
for k = 1%:nfitparams
profindex = pfitID(k); % PROFILE the kth fit parameter 
    
    profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
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
        thetacurr = params(ikeep);
        pboundscurr = pbounds(ikeep);
        pcomb = [pset, currp];
        psetcurr = pcomb(ordp);
        [paramstemp,~, ~, fvaltemp, ~, ~]=fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitIDcurr, psetIDcurr, thetacurr, psetcurr, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pboundscurr);
        
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