% use svd to quickly calculate the Mahalanobis distance
function [t] = quadform(X,sigma)

    len = size(X,1);
    [S, V, D] = svd(sigma);
    newX = X*S;
    newX2 = X*D;
    di = diag(V);
    di = repmat(di',len,1);
    t = sum(newX.*newX2.*di,2);
