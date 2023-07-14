
%% Main file to run on SHARC HPC
clc; clear;

%% Input for each case, change here

% load the historical and sythetic discharge values
O = load('besik_observed.txt');
O_60_50 = load('besik_50x60_years.txt');
O2500 = O_60_50(:,11*365+1:end); %take out the first 10 years
O1000 = reshape(O2500',50*17885,1);
clear O_60_50 O2500

% load the MO solutions as an input to run each design alternatives under uncertain factors
X = load('X200_besik.txt');
X(201,:) = [1.55, 3.5, 3.5, 2, 2];% %identical 10% discount
X(202,:) = [1.55, 3.5, 3.5, 2, 2];% identical


% Project inputs from the DSI Report
HP.hg = 117.3; % initial storage elevation (m), gross head
HP.ht = 0;  % tail water elevation, depth of outflow to stream (m)
HP.L = 208; % the length of penstock (m)
MFD = 0.63; % the minimum environmental flow (m3/s)

HP.cf = 0.15; %the so-called site factor, ranges between 0 and 1.5 (used for the cost of the civil works)
HP.om  =  0.01; % ranges between 0.01 and 0.04,(used for maintenance and operation cost)
HP.fxc  =  5; % the expropriation and other costs including transmission line

%% create SOWs: first future FDCs
S_mod = CalcMetrics(O); % Calculate the summary metrics of the discharge data

[FDC] = CalcFDC(O);

[~,p_pred] = P_rmse(S_mod(:),FDC(:,1),FDC(:,2)); % Now predict function

a = S_mod(1); b = S_mod(2); c = S_mod(3); k = p_pred; clear S_mod

% Derive kp for 1000 years
[B,I] = sort(O1000); %

[FDC2] = CalcFDC(B)  ; % Calculates the Flow Duration Curve
x = FDC2(:,1);

kp =  deal (NaN(numel(x),1));
for m = 1: numel(x)
    if x(m) > c
        kp(m) = erfc(1/(sqrt(2)*b)*log((x(m)-c)/(a-c)))/2;
    else
        kp(m) = 1;
    end
end


%  %return the originial sequence
idx = deal (NaN(numel(O1000),1));
for i = 1:numel(O1000), idx(i,:) = find (I == i);end

O_original = B(idx);
kpn = kp(idx);

%% Load LHS Sampling for all pars

multiplier = load('general_multiplier.txt');

%% Derive FDC a, b and c values based on sampled statistics

% Sampled Parameters
S_median = multiplier (:,5)*median(O1000); % Median
cv  = std(O1000) / mean(O1000);
S_CV = multiplier (:,6)*cv; % Coefficient of variation (CV)
perc1 = prctile(O1000,1); % calculate 1 percntile
S_perc1 = multiplier (:,7)*perc1; % Coefficient of variation (CV)

% check for available futures
if perc1 > 0
 idxm = find(S_median > S_perc1*1.2);
else 
 idxm = find(S_median > S_perc1 + 0.1);
end
    
S_median = S_median(idxm);
S_CV = S_CV(idxm);
S_perc1 = S_perc1(idxm);
multiplier = multiplier(idxm,:);
Ns = size(multiplier,1); %sample of generated data

%% FDC parameter generation
%The climate-perturbed FDC generation model  is accessible from the Zenodo open-access repository at 
%..https://doi.org/10.5281/zenodo.7662679 (Yildiz, 2023), with a link to the 
%...GitHub source codes of the latest release, including a detailed run guide and input files to statistically generate plausible streamflow futures.

Nsize = size(kp,1);

k1st = exp( sqrt(2)*erfcinv(2*0.99));
%initial guess
initial_guess = S_CV;

options = optimset('MaxFunEvals',1000);

pars =  deal (NaN(length(multiplier),1));

for i = 1:Ns
    pars(i) = fminsearch(@(pars) funcFDCgen(pars,kp,S_median(i),S_CV(i),S_perc1(i), Nsize, k1st),initial_guess(i,:),options);
end

bnew = pars;
[cnew, dnew,anew] =  deal (NaN(length(multiplier),1));

for i = 1:Ns
     anew(i) = S_median(i);
     cnew(i) = (S_perc1(i)  - S_median(i)*k1st^bnew(i)) / (1 - k1st^bnew(i));
end

FDCpars = [anew, bnew, cnew ];
%
%% Run the analysis
type = round(X(:,4)); conf = round(X(:,5));

ep = 0.055; % base scenerio
ir = 0.095; % base scenerio
HP.pt = 1500; % steel penstock price per ton ($/kWh)
HP.N = 49;% life time of the project (years)

nf = numel(X(:,1)); %number of alternatives

%%
% constant parameters for all types of projects
HP.maxT = length(kp);           % the size of time steps
HP.e = 0.45*10^(-4);        % epsilon (m) ; fiberglass e = 5*10^(-6) (m), concrete e = 1.8*10^(-4) (m)
HP.v = 1.004*10^(-6);       % the kinematics viscosity of water (m2/s)
HP.g = 9.81;                % acceleration of gravity (m/s2)
HP.ng = 0.98;               % generator-system efficiency
HP.hr = 8760;               % total hours in a year
HP.nf = [0.05 0.33];        % specific spped range of francis turbine
HP.nk = [0.19 1.55];        % specific spped range of kaplan turbine
HP.np = [0.005 0.0612];     % specific spped range of pelton turbine
HP.mf = 0.40;               % min francis turbine design flow rate
HP.mk = 0.20;               % min kaplan turbine design flow rate
HP.mp = 0.11;               % min pelton turbine design flow rate

%% Define variables and interpolation function for calculation of turbines efficiencies

EffCurves = xlsread('EffCurves.xlsx');
HP.perc = EffCurves(:,1);

% Kaplan turbine efficiency
HP.eff_kaplan = EffCurves(:,2);

% Francis turbine efficiency
HP.eff_francis = EffCurves(:,3);

% Pelton turbine efficiency
HP.eff_pelton = EffCurves(:,4);

%%

HP.ny = HP.maxT/365;

[AR, AC, BC, AAE, T_cost, WEP1, power, NPV,PayB] = deal (NaN(Ns,nf));

[Y_ae] = deal(NaN(2450*Ns,nf));

[PB10,  PB15, PB20, nNPV] = deal (NaN(length((kp))/365/50*Ns,nf));

per = 50; %number of realizations


for  j = 1: Ns

    Qsow =  cnew(j) + (anew(j)- cnew(j))*exp(sqrt(2)*bnew(j) *erfcinv(2*kpn));
    Onew = max(Qsow - MFD, 0);
    HP.i = ir*multiplier(j,1);
    HP.ep10 = ep*multiplier(j,2);
    HP.ep40 = ep*multiplier(j,3);
    HP.cor = multiplier(j,4);

    HP.CRF = HP.i*(1+HP.i)^HP.N/((1+HP.i)^HP.N-1); % capital recovery factor
    HP.crf10 = HP.i*(1+HP.i)^10/((1+HP.i)^10-1);
    HP.crf40 = HP.i*(1+HP.i)^39/((1+HP.i)^39-1);
    HP.tf  = 1 / ( 1 + HP.i)^25;

    for k = 1:nf

        if conf(k)==1

            [BC(j,k), power(j,k), AAE(j,k), NPV(j,k), PayB(j,k),  PB15(per*j-(per-1):per*j,k), nNPV(per*j-(per-1):per*j,k)] = calc_hpp1 ( X(k,:) , HP , Onew , type(k),1,k );

        elseif conf(k)==2

            [BC(j,k), power(j,k), AAE(j,k), NPV(j,k), PayB(j,k),  PB15(per*j-(per-1):per*j,k), nNPV(per*j-(per-1):per*j,k)] = calc_hpp2 ( X(k,:) , HP , Onew , type(k), 2, k);

        else

            [BC(j,k), power(j,k), AAE(j,k), NPV(j,k), PayB(j,k),  PB15(per*j-(per-1):per*j,k), nNPV(per*j-(per-1):per*j,k)] = calc_hpp3 ( X(k,:) , HP , Onew , type(k), 3, k);
        end

    end

end

save ('outputBESIK','BC','power', 'AAE', 'NPV', 'PayB','multiplier','FDCpars','kpn','PB15','nNPV','-v7.3')




