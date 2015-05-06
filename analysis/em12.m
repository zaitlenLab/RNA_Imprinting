function [x like] = em12(tmp,x0)

tolx = 1e-3;
toll = 1e-3;

x = x0;
indn = size(tmp,1);
like = sum(logsum_faster(tmp+repmat(log(x),[indn 1])),1);

while (1),
    zs = normrows_faster(tmp+repmat(log(x),[indn 1]));
    prevx = x;
    x = sum(zs,1);
    x = x/sum(x,2);
    prevlike = like;
    like = sum(logsum_faster(tmp+repmat(log(x),[indn 1])),1);
    if (max(abs(x-prevx))<tolx && ((like-prevlike)/prevlike)<toll)
        break;
    end
end
