% calculate sliced Wasserstein distance
function sw = mixturedistancesw(samples, mu, sigma, ppi, ellipse, para_stru, t_vec)   
    [K,M] = size(mu); 
    if (M == 2)
        angle = linspace(0,2*pi,180);
        theta = zeros(2,180);
        theta(1,:) = cos(angle'); theta(2,:) = sin(angle');
    else
        theta = randon_proj(180,M);
    end
    [~,L] = size(theta);
    radon_samples = samples*theta; % projected samples: NxL
    sw = 0;
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
        [~, ~, swt] = optimal_map(pdf_p,pdf_e,(t_vec(end) - t_vec(1))); 
        sw = sw + swt;
    end
    sw = sw/180;