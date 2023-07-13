%The climate-perturbed FDC generation model  is accessible from the Zenodo open-access repository at 
%..https://doi.org/10.5281/zenodo.7662679 (Yildiz, 2023), with a link to the 
%...GitHub source codes of the latest release, including a detailed run guide and input files to statistically generate plausible streamflow futures.
function [Rerr] = funcFDCgen(pars, kp, S_median, S_CV, S_perc1, Nsize, k1st)

b = pars;
p_pred = kp;

d = (S_median - S_perc1) / (1 - k1st^b);
c = S_perc1  - d*k1st^b;

qb = sum(exp(sqrt(2)*b*erfcinv(2*p_pred)))/Nsize;
q2b = sum(exp(sqrt(2)*2*b*erfcinv(2*p_pred)))/Nsize;

Q_CV = sqrt( d^2*(q2b -qb^2) )/(c + d*qb);

% Calculate the OF
Rerr = (Q_CV - S_CV)^2;



