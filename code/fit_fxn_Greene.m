function [pbest,model_N] = fit_fxn_Greene(ydatafit, sigmafit, pfID, psetID, pfit, pset, time, Uvec, lengthvec, N0s) 
%Write a fucntion that fits any combination of the paramters in the
%simplest Greene model with the following parameters:
%p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
% To determine which parameters are set and which ones are fit, pID tells
% us this with a 1 for a parameter that must be fit and a 0 for a parameter
% that is set.

% TEST

% Find pfit
% P = num2cell(pID); 
% Initialize parameters
prop = 0;
rs = 0;
carcap = 0;
alpha = 0;
rr = 0;
ds = 0;
dr = 0;
params = [prop, rs, carcap, alpha, rr, ds, dr];
for i = 1:length(params)
    indset= find(ismember(psetID, i));
    if ~isempty(indset)
    params(i) = pset(indset);
    end
    indfit = find(ismember(pfID,i));
    if ~isempty(indfit)
    params(i) = pfit(indfit);
    end
end
P = num2cell(params);
[prop, rs, carcap, alpha, rr, ds, dr] = deal(P{:});
 
% Define transforms 
% for number of variables in pfit
nfit = length(pfit);
pfxform = @(pval)ones(1,nfit).*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)ones(1,nfit).*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


theta = pfit; 

% write a function that takes set parameters and fit parameters, combines
% them, and runs the model forward for the Uvec and 0s provided
modelfun = @(pfit)simmodelgreene(pfit, time, N0s, pset, Uvec, lengthvec, pfID, psetID); 

loglikelihood = @(phat)sum(log(normpdf(yfxform(ydatafit),yfxform(modelfun(pbxform(phat))), sigmafit)));
  
   
    % Write objective functions for each model
    objfun = @(phat)-loglikelihood(phat); 
    phatbest = fminsearch(objfun, pfxform(theta));
    
    pbest = pbxform(phatbest);

    model_N = modelfun(pbest);
    
    % Consider instead using Levenberg Marquadt algorithm to maximize
    % likelihood
    
    tol = 1e-5;
    % want to minimize your objective function, which in some cases is the
    % sum of your errors, and here is the negative loglikelihood.
    
    %Initialize
    negLL = objfun(pfxform(theta));
end