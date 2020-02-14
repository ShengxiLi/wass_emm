% calculate Euclidean (normal) gradients
function [grad_mu, grad_sigma, grad_ppi] = eucmixgradient(samples, mu, sigma, ppi, ellipse, para_stru, theta, t_vec)
    [K,M] = size(mu); [~,L] = size(theta); grad_ppi = zeros(1,K);
    grad_mu = zeros(K,M); grad_sigma = zeros(M,M,K);
    magnif = (t_vec(end) - t_vec(1))/length(t_vec);
    radon_samples = samples*theta; % projected samples: NxL
    for i = 1:L
        proj_mu = zeros(1,K); proj_sig = zeros(1,K); pdf_p = zeros(1, length(t_vec));
        %Optimise map
        for k = 1:K
            proj_mu(k) = mu(k,:)*theta(:,i);
            proj_sig(k) = theta(:,i)'*sigma(:,:,k)*theta(:,i);   
            if proj_sig(k)<0
                error('Non-positive definite matrix!');
            end
            quadr = (t_vec - proj_mu(k)).^2/proj_sig(k);
            in_log_gener = log_generator(quadr,ellipse(k, :),para_stru{k},1); % M=1 here for projected pdf
            pdf_single = max(exp(in_log_gener - 0.5*log(proj_sig(k))),realmin);
            pdf_single = pdf_single./sum(pdf_single).*ppi(k);
            pdf_p = pdf_p + pdf_single;
        end
        slice_samples = radon_samples(:,i);
        pdf_e = emp_pdf(t_vec, slice_samples);
        pdf_p = pdf_p./sum(pdf_p);
        [~, pote,~] = optimal_map(pdf_p,pdf_e,(t_vec(end) - t_vec(1))); 
        %Optimise parametric distributions
        for k = 1:K
            quadr = (t_vec - proj_mu(k)).^2/proj_sig(k);
            gt = exp(log_generator(quadr, ellipse(k, :), para_stru{k}, 1));
            gtprime = gt.*psi_fun(quadr, ellipse(k, :), para_stru{k}, 1);
            cur_mu = -2*ppi(k)*sum((pote)/sqrt(proj_sig(k)).*gtprime.*((t_vec-proj_mu(k))/proj_sig(k)))*magnif*theta(:,i);
            %cur_mu = cur_mu./sqrt(cur_mu'*cur_mu);
            grad_mu(k,:) = grad_mu(k,:) + cur_mu';
            cur_grad = -ppi(k)*sum((pote)/sqrt(proj_sig(k)^3).*(0.5*gt+gtprime.*((t_vec-proj_mu(k)).^2/proj_sig(k))))*magnif*theta(:,i)*theta(:,i)';  
            %cur_grad = cur_grad./max(1,trace(cur_grad*sigma*cur_grad));   %gradient norm
            grad_sigma(:,:,k) = grad_sigma(:,:,k) + cur_grad;
            grad_ppi(k) = grad_ppi(k) + sum((pote)/sqrt(proj_sig(k)).*gt)*magnif;
        end
    end