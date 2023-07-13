function [RMSE,p_pred] = P_rmse(pars,x,p)
% Now predict function


p_pred4 =  deal (NaN(numel(x),1));
a = pars(1); b = pars(2); c = pars(3);
for m = 1: numel(x)
    if x(m) > c
        p_pred4(m) = erfc(1/(sqrt(2)*b)*log((x(m)-c)/(a-c)))/2;
    else
        p_pred4(m) = 1;
    end
    p_pred = p_pred4;
end


% Calculate RMSE
%RMSE = sqrt ( sum ( ( p_pred - p).^2) / ( numel(x) ) );
RMSE =  sum ( abs( p_pred - p)) / ( numel(x));