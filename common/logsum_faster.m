function [res] = logsum_faster(zs)
% takes a matrix with log probabilities
% returns a vector in which the i'th entry contains the log of the sum of
% probabilities in the i'th row of the input matrix
% (avoids direct exponentiation)

n = size(zs,1);
m = size(zs,2);

tmp = Inf(n,m);
allokay = false;
for k = (1:m),
    t = zs-repmat(zs(:,k),[1 m]);
    tmp(:,k) = zs(:,k)+log(sum(exp(t),2));
    okay = ~(isinf(tmp(:,1:k))|isnan(tmp(:,1:k)));
    if (all(any(okay,2))),
        allokay = true;
        break;
    end
end

if (~allokay),
    disp('error');
    res = [];
    return;
end

[trash ind] = sort(okay,2,'descend');
ind = ind(:,1); % tmp is a vector of class numbers per gene
res = tmp(sub2ind([n m],(1:n)',ind));
