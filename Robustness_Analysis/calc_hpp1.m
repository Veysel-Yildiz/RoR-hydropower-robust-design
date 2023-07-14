function [ BC, power, AAE,  NPV, PayB,  PB15, nNPV] = calc_hpp1 ( x , HP , O  , nt, On ,kr )

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
f_d = moody ( ed , Re_d );

% Claculate flow velocity in the pipe for design head
V_d = 4 * Od / ( pi * D^2 );

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

hf_d = f_d*(HP.L/D)*V_d^2/(2*HP.g)*1.1; % 10% of local losses
%     hl_d = HP.K_sum*V_d^2/(2*HP.g);

HP.hd = HP.hg - hf_d;
power  = HP.hd * HP.g  * Od;
HP.power = power;
HP.qd   = Od;

%%

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

costP = cost_hpp ( HP ,D , Od , 0 ,  nt, On);

% Unpack costs
cost_em = costP(1); cost_pen = costP(2);  cost_ph = costP(4); %tp = costP(3);
%%
if kr == 201 % makes it 201 or 202 
    cost_em = cost_em*0.9; % for identical turbines
end

%% 
cost_cw = HP.cf * (cost_pen + cost_em ); % (in dollars) civil + open channel + Tunnel cost

% Determine total cost (with cavitation)
Cost_other = cost_pen + cost_ph + cost_cw;

T_cost = cost_em * (1+ HP.tf) + (Cost_other + HP.fxc)* HP.cor;

cost_OP = cost_em * HP.om; % operation and maintenance cost

AAE = mean(P) * HP.hr/10^6; % Calculate average annual energy

[Y_p] = deal (NaN(HP.ny,1));
for ny = 1:HP.ny
    Y_p(ny,:) = sum(P(365*ny-364:365*ny));
end

% WEP_1  = prctile(Y_p,1)*24/10^4/AAE; % percentage

R10 = AAE* HP.ep10*0.97; % first 10 years revenue
R40 = AAE* HP.ep40*0.97; % last 40 years revenue

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
    
PB10(j,:) = (sum(Y_p (49*j-9:49*j,:))*24 * 0.97 * HP.ep10/10^6 - cost_OP*10)/Cost10;
      
    PB15(j,:) = ((sum(Y_p (49*j-4:49*j,:)) * HP.ep40 +  sum(Y_p (49*j-14:49*j-5,:)) * HP.ep10) *24 * 0.97/10^6 - cost_OP*15) /Cost15;
    
    PB20(j,:) = ((sum(Y_p (49*j-9:49*j,:)) * HP.ep40  +  sum(Y_p (49*j-19:49*j-10,:)) * HP.ep10) *24 * 0.97/10^6 - cost_OP*20) /Cost20;
    
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



