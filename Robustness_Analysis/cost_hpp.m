function [ cost ] = cost_hpp( HP ,D , Od1 , Od2 , nt, On )


% Thickness of the pipe  [m]
tp  = 8.4/1000*D + 0.002;

cost_pen = pi * tp * D * HP.L * 7.874 * HP.pt/10^6;

% % Calculate the cost of power house (in M dollars)
cost_ph = 200 * (HP.power/1000)^-0.301  * HP.power/10^6;

% Switch among the different turbine combinations

if nt == 2 % Francis turbine cost
    cost_em = 2.927 * (HP.power/1000)^1.174 *(HP.hd)^-0.4933*1.1; % in $
elseif nt == 3 % pelton turbine cost
    cost_em = 1.984 * (HP.power/1000)^1.427 *(HP.hd)^-0.4808*1.1; % in $
else % Kaplan turbine cost
    switch On
        case 1
            cost_em = 2.76 * (HP.power/1000)^0.5774 *(HP.hd)^-0.1193*1.1; % in $
        case 2
            Pn1 =  HP.power* Od1 / HP.qd; Pn2 =  HP.power* Od2 / HP.qd;
            cost_em  = (2.76 * (Pn1/1000)^0.5774 *(HP.hd)^-0.1193*1.15 + 2.76 * (Pn2/1000)^0.5774 *(HP.hd)^-0.1193)*1.1;
        case 3
            Pn1 =  HP.power* Od1 / HP.qd; Pn2 =  HP.power* Od2 / HP.qd;
            cost_em  = (2*2.76 * (Pn1/1000)^0.5774 *(HP.hd)^-0.1193*1.15 + 2.76 * (Pn2/1000)^0.5774 *(HP.hd)^-0.1193)*1.1;
    end
end

cost = [ cost_em , cost_pen, tp, cost_ph ];


