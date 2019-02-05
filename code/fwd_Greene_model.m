function [Nsr, tcrit] = fwd_Greene_model(p, tvec, U, dt, Ncrit)
% This function runs the Greene model of S & R cells with drug-induced
% resistance forward

P = num2cell(p); 
[S0, R0, rs, carcap, alpha, rr, ds, dr] = deal(P{:}); % our parameters
tend = tvec(end);
ttot = 1:dt:tend;
S = zeros([length(ttot),1]);
R = zeros([length(ttot),1]);
N = zeros([length(ttot),1]);
S(1)=S0; 
R(1) = R0; 
N(1) = S(1) + R(1);

for t = 2:tend
    growth_s = rs*S(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    growth_r = rr*R(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    death_s = ds*U(t-1)*S(t-1);
    death_r = dr*U(t-1)*R(t-1);
    trans_s = alpha*U(t-1)*S(t-1);
   
    S(t) = S(t-1) + dt*(growth_s - death_s - trans_s);
    R(t) = R(t-1) + dt *(growth_r - death_r + trans_s);
    N(t) = S(t) + R(t);

end

Nsr = horzcat(N, S, R);
 ikeep = find(ismember(tvec,ttot));
 Nsr= Nsr(ikeep,:);

icrit = find(N>Ncrit,1, 'first');
tcrit= ttot(icrit);

end