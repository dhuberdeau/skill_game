function xx = exclude_outliers(x)

% David Huberdeau, 09/15/2018

xx = x;
temp = nanzscore(x, [], 2) > 3;
xx(temp) = nan;