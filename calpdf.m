function [out_pdf, out_post] = calpdf(in_samples, in_mu, in_sigma, in_pi, ellipse, para_stru)
% out_pdf: output for probability
% out_post: output for posterior
    K=size(in_pi,2);
    [N,M]=size(in_samples);
    in_t = zeros(N,K);
    in_log_gener = zeros(N,K);
    in_detsigma = zeros(1,K);
    in_trsigma = zeros(1,K);
          
    for idx = 1:K
        if (det(in_sigma(:,:,idx))<realmin)
            error('Sigma is singular');
        end             
        in_temp = (in_samples - repmat(in_mu(idx,:), N, 1));
        in_t(:,idx) = max(quadform(in_temp, inv(in_sigma(:,:,idx))),realmin); % NxK
        in_trsigma(idx) = trace(inv(in_sigma(:,:,idx)));
        in_detsigma(idx) = det(in_sigma(:,:,idx))^(0.5);
        in_log_gener(:,idx) = log_generator(in_t(:,idx),ellipse(idx, :),para_stru{idx},M); % NxK
    end
    
    out_pdf = in_log_gener - repmat(log(in_detsigma),N,1); % NxK lisx
    out_pdf = exp(out_pdf);
    out_post = repmat(in_pi,N,1).*out_pdf; 
    location = find(sum(out_post,2)==0); out_post(location,:) = out_post(location,:) + realmin; % plus eps to avoid NAN
    out_post = out_post./repmat(sum(out_post,2),1,K); % \xi NxK
    
    