function [FDC] = CalcFDC(x)
% Calculates the Flow Duration Curve

% Streamflow should be column vector (vertical)
n = size(x,1);

% Now sort the data
x = sort(x);

% And calculate the probabilities
p = ((1:n)-0.5)' ./ n;

% Now return (and we work with exceedance probability), so 1 - p!!
FDC = [x 1 - p];
