% generate random projections
function theta = randon_proj(L,M)
    theta = zeros(M,L);
    thr = 0.99;
    if M==8
        thr = 0.74;
    end
    if M==16
        thr = 0.53;
    end
    for i = 1:L
        ctheta = randn(M,1);
        ctheta = ctheta./sqrt(ctheta'*ctheta);
        theta(:,i) = ctheta;
        if i==1
            continue;
        end
        while max(abs(ctheta'*theta(:,1:(i-1))))>thr
            ctheta = randn(M,1);
            ctheta = ctheta./sqrt(ctheta'*ctheta);
        end
        theta(:,i) = ctheta;
    end