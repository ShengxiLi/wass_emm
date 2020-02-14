% generate synthetic data samples
function [samples, truth_mu, truth_sigma] = systhesis_Gaussiandata(M, K, c, N, e)
samples = [];
    for k = 1:K
        rand_positive = rand(M,M); rand_positive = rand_positive*rand_positive';
        [umatrix, ~] = eig(rand_positive);
        eigval = ones(1,M) + (e - 1).*rand(1,M); eigval(1) = e; eigval(end) = 1; eigval = diag(eigval);
        truth_sigma(:,:,k) = umatrix'*eigval*umatrix;
    end
    distance = (e*c)^0.5;
    for k = 1:K
        truth_mu(k,:) = zeros(1,M) + (k-1)*distance;
        temp = mvnrnd(truth_mu(k,:), truth_sigma(:,:,k), floor(N/K));
        samples = [samples; temp];
    end
    
    
