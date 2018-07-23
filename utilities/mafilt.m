function [fmean, fvar] = mafilt(x, n)
% function f = mafilt(x, n)
% 
% Filters the signal(s) in x by applying a moving average filter with
% non-overlapping bins of size n.  If the length of the signal is not a
% multiple of n, the last bin with partial size is omitted.
%
% INPUT: x - a vector (mx1 or 1xm) or a matrix (m x n-signals) of input
% signals
%
% OUTPUT: fmean - a vector or matrix of length floor(length(x)/n) with
% otherwise same dimensionality as the vector/matrix x with MA values
%         fvar - a vector or matrix of as fmean but with var values. 
%
% Note, get standard error by taking sqrt(fvar)/sqrt(n)
%
% David Huberdeau
% 01/28/14

if min(size(x)) < 2
    % x is vector
    revisedInds = 1:(size(x,1)-mod(size(x,1), n));
    x_mat = reshape(x(revisedInds), n, floor(length(x)/n));
    fmean = mean(x_mat,1);
    fvar = var(x_mat,0,1);
    dimX = size(x);
    if dimX(2) < dimX(1)
        fmean = fmean';
        fvar = fvar';
    end
else
    % x is matrix, dimensions must be set.
    dimX = size(x,2);
    fmean = nan(floor(size(x,1)/n), dimX);
    fvar = nan(floor(size(x,1)/n), dimX);
    revisedInds = 1:(size(x,1)-mod(size(x,1), n));
    for i_col = 1:dimX
        tempx = reshape(x(revisedInds,i_col), n, floor(size(x,1)/n));
        fmean(:,i_col) = mean(tempx,1)';
        fvar(:,i_col) = var(tempx,0,1)';
    end
end






