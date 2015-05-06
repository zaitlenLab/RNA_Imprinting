% target function used by param4model.m to estimate the beta parametesr of
% the imbalanced and imprinted beta distributions
function [tmp] = betamle_23(x,currefs,curalts)

tmp1 = betaln(currefs+x(1),curalts+x(2));
tmp2 = betaln(currefs+x(2),curalts+x(1));
tmpa = tmp1+log(1+exp(tmp2-tmp1)); % discarding /2
tmpb = tmp2+log(1+exp(tmp1-tmp2)); % discarding /2
isnana = isnan(tmpa);
isnanb = isnan(tmpb);
isinfa = isinf(tmpa);
isinfb = isinf(tmpb);
okaya = ~(isnana|isinfa);
okayb = ~(isnanb|isinfb);

if (0<sum(~(okaya|okayb))),
    disp('error');
    tmp = 1
    return;
end
res = tmpa;
res(~okaya) = tmpb(~okaya);

tmp = sum(res)-length(currefs)*betaln(x(1),x(2));
tmp = -tmp;
