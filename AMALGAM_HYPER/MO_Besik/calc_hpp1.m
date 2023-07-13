function [AR, AC,  WEP_1] = calc_hpp1 ( x , HP , O , nt, On  )

% Unpack the parameter values
D = x(1); Od = x(2);

% Pre-allocate output variables
[P] =  deal (NaN(HP.maxT,1));

% Calculate the relative roughness: epsilon / diameter.
ed = HP.e / D;

% design head ------------------------------------------

% Calculate the Reynolds number for design head
Re_d = 4 * Od / ( pi * D * HP.v );

% Find f, the friction factor [-] for design head
f_d = moody (ed , Re_d );

% Claculate flow velocity in the pipe for design head
V_d = 4 * Od / ( pi * D^2 );

if V_d > 9 || V_d < 2.5
    AR = -19999990*V_d ;
    AC = 19999990*V_d;
    WEP_1 = -19999990*V_d; return
end

%%
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

hf_d = f_d*(HP.L/D)*V_d^2/(2*HP.g)*1.1; % 10% of local losses
%     hl_d = HP.K_sum*V_d^2/(2*HP.g);

HP.hd = HP.hg - hf_d;
power  = HP.hd * HP.g  * Od;
HP.power = power;
HP.qd   = Od;


%%
% Now calculate the cavitation costs , suction heads and speeds of turbines

ss_L = 3000/60 * sqrt(Od)/(HP.g*HP.hd)^0.75;
ss_S = 214/60 * sqrt(Od)/(HP.g*HP.hd )^0.75;

if var_name_cavitation(2) <= ss_S  || ss_L <= var_name_cavitation(1)
    
    AR = -29999990 ;
    AC = 29999990;
    WEP_1 = -29999990;return
    
end

perc = HP.perc;
L = HP.L;
ve = HP.v;
hg = HP.hg;
ng = HP.ng;
% % % % % % % % % % % % % Now iterate over time % % % % % % % % % % % % % %
parfor t = 1 : HP.maxT
    
    % Check sum of Od1 and Od2
    q = min (O(t) , Od);
    
    if  q < kmin*Od
        % Both turbines shut down
        P(t) = 0;
        
    else
        % Calculate the Reynolds number
        Re = 4 * q / ( pi * D * ve );
        
        % Find f, the friction factor [-]
        f  = moody ( ed , Re );
        
        % Claculate flow velocity in the pipe
        V = 4 * q / ( pi * D^2 );
        
        % Calculate the head loss due to friction in the penstock
        hf = f *(L/D)*V^2/(2*9.81)*1.1; % hf=f*(HP.L/D)*(V^2/2g),
        
        hnet = hg - hf;
        
        % large francis/kaplan/cross turbine efficiency
        %         n1 =  eff( q , Od ) ;
        ck = q/Od;
        n = interp1(perc,func_Eff,ck);
        
        P(t) = hnet * q * 9.81 * n * ng;
        
    end
end

% % % % % % % % % % % % % End iterate over time % % % % % % % % % % % % % %

% Determine the total costs and profit of operation
costP = cost_hpp ( HP ,D , Od , 0 , nt, On );

% Unpack costs
cost_em = costP(1); cost_pen = costP(2);  cost_ph = costP(4); %tp = costP(3);
cost_cw = HP.cf * (cost_pen + cost_em ); % (in dollars) civil + open channel + Tunnel cost

%%
% Determine total cost (with cavitation)
Cost_other = cost_pen + cost_ph + cost_cw;

T_cost = cost_em * (1+ HP.tf) + Cost_other + HP.fxc;

cost_OP = cost_em * HP.om; % operation and maintenance cost

AAE = mean(P) * HP.hr/10^6; % Calculate average annual energy

[Y_p] = deal (NaN(HP.ny,1));
for ny = 1:HP.ny
    Y_p(ny,:) = sum(P(365*ny-364:365*ny));
end

WEP_1  = prctile(Y_p,1)*24/10^6; % percentage

AR = AAE* HP.ep*0.97; % AnualRevenue in M dollars 5% will not be sold

AC = HP.CRF * T_cost + cost_OP; % Anual cost in M dollars

if AR/ AC  < 0.90
    AC =  AC*10; AR = AR/10; WEP_1= WEP_1/10;
end


