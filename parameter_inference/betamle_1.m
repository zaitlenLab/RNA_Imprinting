% target function used by param4model.m to estimate the beta parametesr of
% the balanced beta distribution
function [tmp] = betamle_1(x,currefs,curalts)

tmp = betaln(currefs+x(1),curalts+x(2));
tmp = sum(tmp)-length(currefs)*betaln(x(1),x(2));
tmp = -tmp;
