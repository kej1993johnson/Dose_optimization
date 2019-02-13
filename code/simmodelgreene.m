function [Y] = simmodelgreene(pfit,ytimefit, N0s, pset, Uvec, lengthvec)
% SIMMODELGreene returns Y values corresponding to times T for parameter vector p
% for James Greene model
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R ;
% 
%pset = [S0, R0, rs, carcap];
%pfit = [ alpha, rr, ds, dr];

% define inputs
P = num2cell(pfit); 
[alpha, rr, ds] = deal(P{:}); % our parameters

P2 = num2cell(pset);
[rs, carcap, props, dr] = deal(P2{:});




Nsr = [];
tsum = [];
tdrug = 1; % since we are monitoring first treatment only
istart = vertcat(1,cumsum(lengthvec(1:end-1,1))+1);
iend = cumsum(lengthvec(:,1));

ist = vertcat(1,cumsum(lengthvec(1:end-1,2))+1);
ie = cumsum(lengthvec(:,2));

for i = 1:size(lengthvec,1)
    S0 = props*N0s(i);
    R0 = (1-props)*N0s(i);
    % vary these based on what we're fitting
  
    p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
    
    tvec = round(ytimefit(istart(i):iend(i)),0);
    U = Uvec(ist(i):ie(i));
    dt = 1;
    [Nsri, ~, ~] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    
    Nsr = vertcat(Nsr, Nsri);
end



Y = Nsr(:,1);
end