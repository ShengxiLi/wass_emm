% generate synthetic data samples with all random cases
function [samples, truth_mu, truth_sigma, truth_ppi] = systhesis_Gaussiandata_all_random(M, K, c, N, e, input_rand_sigma)
    samples = [];
    pi_list = rand(1,K);
    pi_list = pi_list./(sum(pi_list))*N;
    pi_list = floor(pi_list);
    pi_list(end) = N - sum(pi_list(1:end-1));
    truth_ppi = pi_list./N;
    assert(~(e < 1), 'e should be larger than 1 cause e = \lambda_{max}/lambda_{min}.');
    for k = 1:K
       rand_positive = rand(M,M); rand_positive = rand_positive*rand_positive';
       [umatrix, ~] = eig(rand_positive);
%         [umatrix, ~] = eig(input_rand_sigma(:,:,k));
        eigval = ones(1,M) + (e - 1).*rand(1,M); eigval(1) = e; eigval(end) = 1; eigval = diag(eigval);
        eigvalue = diag(eigval); eigvalue = eigvalue./mean(eigvalue);
        eigval = diag(eigvalue);
        truth_sigma(:,:,k) = umatrix'*eigval*umatrix;
    end
    setting_distance = e*c;    
    truth_mu(1,:) = (rand(1,M)-0.5) * 200;
    temp = mvnrnd(truth_mu(1,:), truth_sigma(:,:,1), pi_list(1));
    samples = [samples; temp];
    for k = 2:K
        min_distance = inf;
        while ((min_distance < setting_distance)||(min_distance > 2*setting_distance))
            min_distance = inf;
            for search_idx = 1:(k-1)
                choosing_dir = rand(1,M) - 0.5;
                choosing_dir = choosing_dir./(norm(choosing_dir)).*setting_distance*(rand(1,1)+1);
                rp = randperm(k-1);
                truth_mu(k,:) = truth_mu(rp(1),:) + choosing_dir;
                curr_distance = norm(truth_mu(k,:) - truth_mu(search_idx,:));
                if min_distance > curr_distance
                    min_distance = curr_distance;
                end
            end
        end
        temp = mvnrnd(truth_mu(k,:), truth_sigma(:,:,k), pi_list(k));
        samples = [samples; temp];
    end