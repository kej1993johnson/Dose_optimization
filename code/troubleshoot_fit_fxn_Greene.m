%[pbest,model_N] = fit_fxn_Greene(ydatafit, sigmafit, pfID, psetID, pfit, pset, time, Uvec, lengthvec, N0s)
%ydatafit,sigmafit, pfID, psetID, pfitguess, pset, ytimefit, Uvec, lengthvec, N0s
pfit = pfitguess;
time = ytimefit;

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
    
    %% Consider instead using Levenberg Marquadt algorithm to maximize
    % likelihood
    
    tol = 1e-5;
    lambda =2;
    it_count = 0;
    negLL_c = 0;
    % want to minimize your objective function, which in some cases is the
    % sum of your errors, and here is the negative loglikelihood.
    max_iters = 1;
    deltap = 1e-8;
    spant = length(time);
    %Initialize best guess
    negLLbest = objfun(pfxform(theta));
    negLLinit = negLLbest;
    Nbest = modelfun(theta);
    Ntinit = modelfun(theta);
    pe = theta;
    pbounds = [0, 1; 0, 1; 0,1]; % set your parameter bounds pretty widely
    pevec = [];
    delvec = [];
    weightvec = 1./(sigmafit.^2);
    W =diag(weightvec);
    eta = 1e-5;
    %%
    while negLLbest>tol && it_count<max_iters
        it_count = it_count +1;
        
        % Calculate Jacobian
        for n = 1:nfit
            ptemp = pe;
            ptemp(n) = ptemp(n) + deltap;
            Nt = modelfun(ptemp); % find temporary model output
            
            J(1:spant, n) = (Nt-Nbest)/deltap;% this is the derivative of the model function wrt the parameters
            % This subtracts the current parameter prediction of the N(t)
            % with the best parameter prediction of N(t) at each time
            % point. We only output time points that we have data for, so
            % this doesn't need to be written in terms of the discrete time
            % points we want from the model prediction
        end
            % Calculate change in parameters
            err_v =(ydatafit - Nbest); % evaluate objective function from current best parameters
            negLL_v = objfun(pfxform(pe));
            hgd = eta*J'*W*err_v;
          
            del = (J'*W*err_v)\(J'*W*J + lambda*(diag(J'*W*J)));
            %del = (J'*J+lambda*diag(diag(J'*J)))\((J'*err_v));
            delvec = vertcat(delvec, del');
            
            % Update parameters by adding del
            % reset phattemp to pe (back to your current parameters that
            % were probed while calculating the Jacobian
            ptemp = pe;
            pevec= vertcat(pevec,pe);
            % transform into reals to impose real bounds 
            for n = 1:nfit
                ptemp(n)= ptemp(n)+del(n);
                if ptemp(n)< pbounds(n,1)
                    ptemp(n) = pbounds(n,1);
                elseif ptemp(n)>pbounds(n,2)
                    ptemp(n) = pbounds(n,2);
                end
            end
            
            %Now have an updated phattemp
            % Evaluate the model
            Nt = modelfun(ptemp);
            negLLt = objfun(pfxform(ptemp));
            if negLLt < negLLbest
                negLL_c(it_count) = negLLt;
                negLLbest = negLLt;
                pe = ptemp;
                lambda = lambda/2;
                eta = eta/2;
                Nbest = Nt;
            else
                negLL_c(it_count) = negLLt;
                lambda = lambda*4;
                eta = eta*4;
            end
    end
    %%
    
    figure;
    imagesc(J)
    
    pbestLM = pe
 figure;
 plot(time, Ntinit, 'g*', 'LineWidth',3)
 hold on
 plot(time, Nt, 'b*', 'LineWidth',1)
 hold on
 plot(time, Nbest, 'r.', 'LineWidth',2)
 plot(time, ydatafit, 'k*', 'LineWidth', 2)
 xlabel('time (hours)')
 ylabel('N(t)')
 legend('Initial Guess', 'current N(t)', 'best N(t)', 'real data')
 
 figure;
 plot(1:1:max_iters, negLL_c, 'b*')
 hold on
 plot(1, negLLinit, 'r*')
 legend('current negLL', 'initial negLL')
 xlabel('iteration')
 ylabel('objective function (negative log likelihood)')
 
 figure;
 subplot(1,2,1)
 plot(time, err_v, '*')
 subplot(1,2,2)
 plot(time, ydatafit, '*')
 hold on
 plot(time, Nbest,'.')
 legend('data', 'current model')
 

 
 
 