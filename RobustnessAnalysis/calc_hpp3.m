function [Revenue, Cost, WEP_1, BC, power, AAE,  T_cost, NPV, PayB, Y_p, PB10,  PB15, PB20, nNPV] = calc_hpp3 ( x , HP , O , nt, On,k )

% Unpack the parameter values

D = x(1); Od1 = x(2); Od2 = x(3);

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
%%



if nt == 2 % Francis turbine
    kmin = HP.mf;
    func_Eff = HP.eff_francis;
    
elseif nt == 3 % Pelton turbine
    kmin = HP.mp;
    func_Eff = HP.eff_pelton;
else
    kmin = HP.mk;
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
costP = cost_hpp ( HP ,D , Od1 , Od2 ,  nt, On);
% Unpack costs
cost_em = costP(1); cost_pen = costP(2);  cost_ph = costP(4); %tp = costP(3);

%%
if k == 201 % makes it 201 or 202 
    cost_em = cost_em*0.9; % for identical turbines
end

%% 
cost_cw = HP.cf * (cost_pen + cost_em ); % (in M dollars) civil + open channel + Tunnel cost

Cost_other = cost_pen + cost_ph + cost_cw;

T_cost = cost_em * (1+ HP.tf) + (Cost_other + HP.fxc) * HP.cor;

cost_OP = cost_em * HP.om ; % operation and maintenance cost

AAE = mean(P) * HP.hr/10^6; % Calculate average annual energy

% WEP_1 = min(WEP)* HP.hr/10^6;
[Y_p] = deal (NaN(HP.ny,1));
for ny = 1:HP.ny
    Y_p(ny,:) = sum(P(365*ny-364:365*ny));
end

WEP_1  = prctile(Y_p,1)*24/10^4/AAE; % percentage

R10 = AAE* HP.ep10*0.98; % first 10 years revenue
R40 = AAE* HP.ep40*0.98; % last 40 years revenue

AR10 = R10/HP.crf10;
AR40 = R40/HP.crf40 * (1+HP.i)^-10;

Revenue = AR10 + AR40;

pb10 = T_cost / (R10 - cost_OP);
if pb10 < 10.01
PayB = pb10;
else 
PayB = (T_cost - 10*(R10 - cost_OP)) / (R40 - cost_OP) + 10;
end

Cost = T_cost + cost_OP / HP.CRF; % Anual cost

BC = Revenue / Cost;

NPV = Revenue - Cost;

Cost10 = cost_em + (Cost_other + HP.fxc) * HP.cor ; %10 year cost
Cost15 = cost_em + (Cost_other + HP.fxc) * HP.cor ; %20 year cost
Cost20 = cost_em + (Cost_other + HP.fxc) * HP.cor ; %20 year cost


[PB10,  PB15, PB20, nNPV] = deal (NaN(50,1));


for j = 1:50
    
    PB10(j,:) = (sum(Y_p (49*j-9:49*j,:))*24 * 0.98 * HP.ep10/10^6 - cost_OP*10)/Cost10;
      
    PB15(j,:) = ((sum(Y_p (49*j-4:49*j,:)) * HP.ep40 +  sum(Y_p (49*j-14:49*j-5,:)) * HP.ep10) *24 * 0.98/10^6 - cost_OP*15) /Cost15;
    
    PB20(j,:) = ((sum(Y_p (49*j-9:49*j,:)) * HP.ep40  +  sum(Y_p (49*j-19:49*j-10,:)) * HP.ep10) *24 * 0.98/10^6 - cost_OP*20) /Cost20;
    
end

for j = 1:50
    
    Rev10 = 0;
    for k =  1:10 %
        Rev10 =  Rev10 + sum(Y_p(49*j-49 + k)*(1+HP.i)^-k) * 24 * 0.98 * HP.ep10/10^6; %10 year doscounted rev
    end
    
    Rev40 = 0;
    for k = 11:49
        Rev40 = Rev40 + sum(Y_p(49*j-49 + k)*(1+HP.i)^-k) * 24 * 0.98 * HP.ep40/10^6;%39 year doscounted rev
    end
    
    nNPV(j,:) = Rev10 + Rev40 - T_cost;
    
end


