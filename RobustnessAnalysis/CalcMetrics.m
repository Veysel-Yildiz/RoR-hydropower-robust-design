function [S_mod] = CalcMetrics(Y)
% Calculate the summary metrics of the discharge data

S_mod = deal (NaN(3,1));

% Derive FDC
FDC = CalcFDC(Y);

% The second column is the exceedance probability
p = FDC(:,2);

% Now define the streamflow
x = FDC(:,1);

options = optimset('MaxFunEvals',1000);

pars = [0.5 1.5 0.01]; %
pars = fminsearch(@(pars) P_rmse(pars,x,p),pars,options);
S_mod(1) = pars(1); S_mod(2) = pars(2); S_mod(3) = pars(3);

