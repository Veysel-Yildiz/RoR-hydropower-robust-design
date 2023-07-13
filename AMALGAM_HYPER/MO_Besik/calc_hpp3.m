function [AR,  AC,  WEP_1] = calc_hpp3 ( x , HP , O , nt, On )

% Unpack the parameter values
D = x(1); Od1 = x(2); Od2 = x(3);

if Od1 < Od2
    AR = -9999990*(Od2-Od1);
    AC = 9999990*(Od2-Od1);
    WEP_1 = -9999990*(Od2-Od1); return
end

% Pre-allocate output variables
[P] =  deal (NaN(HP.maxT,1));

% Calculate the relative roughness: epsilon / diameter.
ed = HP.e / D;

% design head ------------------------------------------

% Calculate the Reynolds number for design head
Re_d = 4 * (2*Od1 + Od2) / ( pi * D * HP.v );

% Find f, the friction factor [-] for design head
f_d = moody ( ed , Re_d );

% Claculate flow velocity in the pipe for design head
V_d = 4 * (2*Od1 + Od2)  / ( pi * D^2 );

if V_d > 9 || V_d < 2.5
    AR = -19999990*V_d ;
    AC = 19999990*V_d;
    WEP_1 = -19999990*V_d; return
end



if nt == 2 % Francis turbine
    kmin = HP.mf;
    var_name_cavitation = HP.nf; %specific speed range
    func_Eff = HP.eff_francis;
    
elseif nt == 3 % Pelton turbine
    kmin = HP.mp;
    var_name_cavitation = HP.np; %specific speed range
    func_Eff = HP.eff_pelton;
else
    kmin = HP.mk;
    var_name_cavitation = HP.nk; %specific speed range
    func_Eff = HP.eff_kaplan;
end

%%

% head losses
hf_d = f_d*(HP.L/D)*V_d^2/(2*HP.g)*1.1;

HP.hd = HP.hg  - hf_d;

qd   = 2*Od1  + Od2; HP.qd   = qd ;
power =  HP.hd * HP.g * qd;
HP.power = power;

%%
% Now calculate the cavitation costs , suction heads and speeds of turbines

ss_L1 = 3000/60 * sqrt(Od1)/(HP.g*HP.hd)^0.75;
ss_S1 = 214/60 * sqrt(Od1)/(HP.g*HP.hd )^0.75;

ss_L2 = 3000/60 * sqrt(Od2)/(HP.g*HP.hd)^0.75;
ss_S2 = 214/60 * sqrt(Od2)/(HP.g*HP.hd )^0.75;

if var_name_cavitation(2) <= ss_S1  || ss_L1 <= var_name_cavitation(1)
    
    AR = -29999990 ;
    AC = 29999990;
    WEP_1 = -29999990;return
end

if var_name_cavitation(2) <= ss_S2  || ss_L2 <= var_name_cavitation(1)
    
    AR = -39999990 ;
    AC = 39999990;
    WEP_1 = -39999990; return
end

perc = HP.perc;
L = HP.L;
ve = HP.v;
hg = HP.hg;
ng = HP.ng;

% % % % % % % % % % % % % Now iterate over time % % % % % % % % % % % % % %
parfor t = 1 : HP.maxT
    
    % Check sum of Od1 and Od2
    q = min (O(t) , qd);
    
    % Calculate the Reynolds number
    Re = 4 * q / ( pi * D * ve);
    
    % Find f, the friction factor [-]
    f  = moody (  ed , Re );
    
    % Claculate flow velocity in the pipe
    V = 4 * q / ( pi * D^2 );
    
    % Calculate the head loss due to friction in the penstock
    hf = f *(L/D)*V^2/(2*9.81)*1.1; % hf=f*(HP.L/D)*(V^2/2g),
    %hf(t,:)= 10.29*0.012^2*O(t)^2*HP.L/(D^5.333); % fiction losses based on manning equation
    
    hnet = hg  - hf;
    
    % Now calculate power (kW) for different turbine combinations
    
    % Now check whether turbines shut down
    if q < kmin*Od2  % O(t) < HP.m2*Od2
        % Both turbines shut down 
        P(t) = 0;
        
    elseif q >= kmin * Od2 && q <=  Od2 % only the small turbine in operation
        
        % small francis/kaplan/cross turbine efficiency
        ck = q/Od2;
        n2 = interp1(perc,func_Eff,ck);
        P(t) = hnet * q * 9.81 * n2 * ng;
        
    elseif q > Od2 && q < Od1 + kmin * Od2 %:  only one turbine in operation, whihcever achives best production
        
        % large francis/kaplan/cross turbine efficiency
        
        q1 = min(Od1, q);
        ck1 = q1/Od1;
        n1 = interp1(perc,func_Eff,ck1);
        P1 = hnet * q1 * 9.81 * n1 * ng;
        
        % small turbine at full capacity.
        n2 = func_Eff(end);
        P2 = hnet * Od2 * 9.81 * n2 * ng;
        % [kW] maximum power produced, monthly
        P(t) =  max( P1, P2 );
        
        %elseif  q >  Od1 + kmin * Od2 && q <= Od1 + Od2 % both turbines in operation eval(evalstr.condition_3)
    elseif q >=  Od1 + kmin * Od2 && q < Od1 + Od2 + kmin * Od1
        % large francis/kaplan/cross turbine efficiency
        n1 =  func_Eff(end);
        P1 = hnet * Od1 * 9.81 * n1 * ng;
        
        % Check flow
        q2 = min(Od2, (q- Od1));
        ck2 = q2/Od2;
        n2 = interp1(perc,func_Eff,ck2);
        P2 = hnet * q2 * 9.81 * n2 * ng;
        
        P(t) = P1 + P2;
        
    elseif q >=  Od1 + Od2 + kmin * Od1 && q < 2*Od1 + kmin * Od2 % three turbines in opoeration
        
        n1 =  func_Eff(end);
        P1 = hnet * Od1 * 9.81 * n1 * ng;
        
        P2 = hnet * Od2 * 9.81 * n1 * ng;
        
        q3 = min(Od1, (q - Od1-Od2));
        ck3 = q3/Od1;
        n3 = interp1(perc,func_Eff,ck3);
        P3 = hnet * q3 * 9.81 * n3 * ng;
        
        P(t) = P1 + P2 + P3;
        
    elseif q >=  2*Od1 + kmin * Od2 %&& q <= 2*Od1 + Od2
        
        n1 =  func_Eff(end);
        P1 = hnet * Od1 * 9.81 * n1 * ng;
        
        q2 = min(Od2, (q-2*Od1));
        
        ck2 = q2/Od2;
        n2 = interp1(perc,func_Eff,ck2);
        P2 = hnet * q2 * 9.81 * n2 * ng;
        
        P(t) = 2*P1 + P2;
        
    end
end


% % % % % % % % % % % % % End iterate over time % % % % % % % % % % % % % %
% Determine the total costs and profit of operation
costP = cost_hpp ( HP ,D , Od1 , Od2 ,  nt, On );
% Unpack costs
cost_em = costP(1); cost_pen = costP(2);  cost_ph = costP(4); %tp = costP(3);
cost_cw = HP.cf * (cost_pen + cost_em ); % (in M dollars) civil + open channel + Tunnel cost

%%
% Determine total cost (with cavitation)
Cost_other = cost_pen + cost_ph + cost_cw;
T_cost = cost_em * (1+ HP.tf) + Cost_other + HP.fxc;

cost_OP = cost_em * HP.om ; % operation and maintenance cost
%%

AAE = mean(P) * HP.hr/10^6; % Calculate average annual energy

[Y_p] = deal (NaN(HP.ny,1));
for ny = 1:HP.ny
    Y_p(ny,:) = sum(P(365*ny-364:365*ny));
end

WEP_1  = prctile(Y_p,1)*24/10^6; % percentage

AR = AAE* HP.ep*0.98; % AnualRevenue in M dollars 5% will not be sold

AC = HP.CRF * T_cost + cost_OP; % Anual cost in M dollars

%BC = AR / AC;
if AR/ AC  < 0.90
    AC =  AC*10; AR = AR/10; WEP_1= WEP_1/10;
end

