function [zs] = normrows_faster(zs)
% takes a matrix with log probabilities
% returns a probabilities matrix in which each rows sums to 1

n = size(zs,1);
m = size(zs,2);

tmp = logsum_faster(zs);
zs = zs-repmat(tmp,[1 m]);
zs = exp(zs);
