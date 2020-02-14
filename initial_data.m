% initialisation
function [ini_mu, ini_sigma, ini_pi] = initial_data(samples, K, Is_kmeans)
    [N,M] = size(samples);
    ini_sigma = zeros(M,M,K);
    ini_pi = zeros(1,K);
    if Is_kmeans
        % Use k-means to initialise, you may need the vl_feat toolkit
        X = samples';
        [mu, ~] = vl_kmeans(X, K, ...
        'Initialization', 'plusplus');%, ...
    %    'MaxNumIterations',5);  
        ini_mu = mu';
    else
        r = randperm(N);
        ini_mu = samples(r(1:K),:) + randn(K,M);
    end

    for i = 1:K
        ini_pi(i) = 1/K;
        ini_sigma(:,:,i) = eye(M); 
    end