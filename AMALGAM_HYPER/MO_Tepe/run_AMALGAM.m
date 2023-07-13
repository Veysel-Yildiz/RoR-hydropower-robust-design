

%% Define inputs, % Model inputs, project-based parameters  that change from one another

clc;
clear;

%% Model inputs, project-based parameters  that change from one another

Q = load('tepe1000years.txt');
%  HP.MFD = 0; % already considered 0.662; % the minimum environmental flow (m3/s)
%

HP.hg = 56; % gross head
HP.ht = 0;
HP.L = 54; % the length of penstock (m)

HP.cf = 0.15; %the so-called site factor, ranges between 0 and 1.5 (used for the cost of the civil works)
HP.om  =  0.01; % ranges between 0.01 and 0.04,(used for maintenance and operation cost)
HP.fxc  =  2; % the expropriation and other costs including transmission line

HP.ep = 0.055; % electricity price in Turkey ($/kWh)
HP.pt = 1500; % steel penstock price per ton ($/kWh)
HP.i = 0.095;% the investment discount rate (or interest rate, %)
HP.N = 49;% life time of the project (years)

%%
HP.maxT = length(Q(:,1));           % the size of time steps
HP.e = 0.45*10^(-4);         % epsilon (m) ; fiberglass e = 5*10^(-6) (m), concrete e = 1.8*10^(-4) (m)
HP.v = 1.004*10^(-6);       % the kinematics viscosity of water (m2/s)
HP.g = 9.81;                % acceleration of gravity (m/s2)
HP.ng = 0.98;                % generator-system efficiency
HP.hr = 8760;               % total hours in a year
HP.nf = [0.05 0.33];        % specific spped range of francis turbine
HP.nk = [0.19 1.55];        % specific spped range of kaplan turbine
HP.np = [0.005 0.0612];       % specific spped range of pelton turbine
HP.mf = 0.40;               % min francis turbine design flow rate
HP.mk = 0.20;               % min kaplan turbine design flow rate
HP.mp = 0.11;               % min pelton turbine design flow rate

%% Define variables and interpolation function for calculation of turbines efficiencies

EffCurves = readmatrix('EffCurves.xlsx');
HP.perc = EffCurves(:,1);
% Kaplan turbine efficiency
HP.eff_kaplan = EffCurves(:,2);

% Francis turbine efficiency
HP.eff_francis = EffCurves(:,3);

% Pelton turbine efficiency
HP.eff_pelton = EffCurves(:,4);

%%

HP.CRF = HP.i*(1+HP.i)^HP.N/((1+HP.i)^HP.N-1); % capital recovery factor
HP.tf  = 1 / ( 1 + HP.i)^25;

HP.ny = HP.maxT/365;
%HP.ny = max(HP.ny,1);
WEP =  deal (NaN(max(HP.ny,1),1));

%% Define fields of AMALGAMPar
AMALGAMPar.N = 100;                         % Define population size
AMALGAMPar.T = 1000;                          % How many generations?
AMALGAMPar.d = 5;                          % How many parameters?
AMALGAMPar.m = 3;                           % How many objective functions?

%% Define fields of Par_info
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'bound';           % Explicit boundary handling

[ti, tj] = select_turbine (HP);

qt = max(min(Q),0.5);
dmin = sqrt(qt/(9*3.14)*4);

% Define parameter ranges, if 'latin'
Par_info.min = [dmin qt 0.3 ti 0.501];  % If 'latin', min values
Par_info.max = [5.00 20 5.0 tj 3.499];   % If 'latin', max values

func_in.Q = Q;
func_in.HP = HP;

% Define name of function
Func_name = 'MODEL_main';
%% Define structure options
options.print = 'no';                      % Print output to screen (tables and figures)
options.parallel = 'yes';
%% Run the AMALGAM code and obtain non-dominated solution set
[X,F,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options,func_in);
