% The core part to solve EMMs via manifold optimisation on an approximate
% Wasserstein distance

% Input: samples: input data (n x m)
%        ini_mu: initial mu (k x m)
%        ini_sigma: initial sigma (m x m x k)
%        ini_pi: initial pi (1 x k)
%        ellipse: specify the type of EMM
%        para_stru: give specific parameters for each mixture
%        opt: optional parameters 

% Output: final_mu: optimised mu
%         final_sigma: optimised sigma
%         final_pi: optimised pi
%         info: status recorded for each check period
%         labels: cluster labels
%         samples_rec: cluster reconstruction by mu of each cluster

function [final_mu, final_sigma, final_pi, info, labels, samples_rec] = mixture_solver(samples, ini_mu, ini_sigma, ini_pi, ellipse, para_stru, opt)   
    K=size(ini_pi,2); [N,M]=size(samples); ite = 0; 
    mu = ini_mu; sigma = ini_sigma; ppi = ini_pi;
  
%% options setting
    % verbosity
    if isfield(opt, 'Is_visulisation')
        Is_visulisation = opt.Is_visulisation;
    else
        Is_visulisation = 2;
    end

    % maximal iterations
    if isfield(opt, 'maxite')
        maxite = opt.maxite;
    else
        maxite = 1000;
    end

    % check period
    if isfield(opt, 'cp')
        checkperiod = opt.cp;
    else
        checkperiod = 1;
    end   

    % learning rate
    if isfield(opt, 'lr')
        lr = opt.lr;
    else
        lr = 0.001;
    end 
    
    % number of random projections per epoch   
    if isfield(opt, 'L')
        L = opt.L;
    else
        L = 1;
    end 
    
    % number of iterations after which to update random projection
    if isfield(opt, 'repnum')
        repnum = opt.repnum;
    else
        repnum = 1;
    end 
    
    % precision
    if isfield(opt, 't_vec')
        t_vec = opt.t_vec;
    else
        emperi_max = max(max(samples))*sqrt(2*M);
        intvals = 1e4; t_vec = linspace(-emperi_max,emperi_max,intvals);
    end 

%% main part for manifold optimisation in a stochastic manner 
    info.time = [];
    info.LL = zeros(1,maxite);   
    info.SW = zeros(1,maxite); 
    weights_man = sqrt(ppi);
    sw = mixturedistancesw(samples, mu, sigma, ppi, ellipse, para_stru, t_vec);
    info.SW(1) = sw;
    [pdf, ~] = calpdf(samples, mu, sigma, ppi, ellipse, para_stru); 
    ll=-mean(log(ppi*pdf'));
    info.LL(1) = ll;
    info.time = [info.time 0];
    theta = randon_proj(L,M); 
    fst_grad_sig = zeros(M,M,K);
    sec_grad_sig = zeros(M,M,K);
    adlr_pro = zeros(1,K);       
    beta1_sig = 0.9;
    beta2_sig = 0.999;

    while ite<maxite    
        ite = ite + 1;
        info.mu{ite} = mu;
        info.sigma{ite} = sigma;
        info.ppi{ite} = ppi;
        info.theta{ite} = theta;
        start_time = tic();
        [grad_mu, grad_sigma, grad_ppi] = eucmixgradient(samples, mu, sigma, ppi, ellipse, para_stru, theta, t_vec);
        grad_ppi = grad_ppi/max(max(abs(grad_ppi)),1);
        grad_ppi = grad_ppi*2.*weights_man;
        rgrad_ppi = grad_ppi - (weights_man(:)'*grad_ppi(:))*weights_man;
        step_ppi = -0.1*lr*rgrad_ppi;
        weights_man = (weights_man + step_ppi)/norm(weights_man + step_ppi);
        ppi = weights_man.^2;  
        for k = 1:K 
            mu(k,:) = mu(k,:) - 10*lr.*grad_mu(k,:)/L; 
            % Dadam part
            fst_grad_sig(:,:,k) = beta1_sig*fst_grad_sig(:,:,k) + (1-beta1_sig)*grad_sigma(:,:,k);
            rgrad = grad_sigma(:,:,k); %*sigma(:,:,k) + sigma(:,:,k)*grad_sigma(:,:,k);
            sec_grad_sig(:,:,k) = beta2_sig*sec_grad_sig(:,:,k) + (1-beta2_sig)*(rgrad*rgrad');
            adlr = max(theta'*sec_grad_sig(:,:,k)*theta,adlr_pro(k));
            adlr_pro(k) = adlr;
            lr_adp = lr*sqrt(1-beta2_sig^ite)/(1-beta1_sig^ite);
            stepsig = -lr_adp*fst_grad_sig(:,:,k)/(sqrt(adlr_pro(k))+eps);
            sigma(:,:,k) = (stepsig+eye(M))*sigma(:,:,k)*(stepsig+eye(M));
        end 

        timecost = info.time(end) + toc(start_time);
        info.time = [info.time timecost];
        if ((mod(ite, checkperiod) == 0)||(ite == maxite))
            if (Is_visulisation == 2)
                sw = mixturedistancesw(samples, mu, sigma, ppi, ellipse, para_stru, t_vec);
                info.SW(ite) = sw;
            end
            [pdf, ~] = calpdf(samples, mu, sigma, ppi, ellipse, para_stru); 
            ll=-mean(log(ppi*pdf'));
            info.LL(ite) = ll; 
            if ~(isreal(ll))
                error('Likelihood is not a real number');
            end  
            if (Is_visulisation == 2)
                disp([num2str(ite) ': Likelihood ' num2str(ll) ' SW ' num2str(sw)]);
            elseif (Is_visulisation == 1)
                disp([num2str(ite) ': Likelihood ' num2str(ll)]);
            end
        end     
        if mod(ite, repnum) == 0
            theta = randon_proj(L,M);
        end        
    end  
    % Do interplotation to plot a smooth curve
    temp_time = info.time;
    info.time = temp_time(2:end);       
    true_value_locate = find(info.LL~=0);
    for i = 1:(length(true_value_locate)-1)
        interval = true_value_locate(i+1) - true_value_locate(i);
        inter_x = 1:(interval+1);
        inter_approx = interp1([1 (interval+1)], [info.LL(true_value_locate(i)) info.LL(true_value_locate(i+1))],inter_x);
        info.LL(true_value_locate(i):true_value_locate(i+1)) = inter_approx;
    end 
    true_value_locate = find(info.SW~=0);
    for i = 1:(length(true_value_locate)-1)
        interval = true_value_locate(i+1) - true_value_locate(i);
        inter_x = 1:(interval+1);
        inter_approx = interp1([1 (interval+1)], [info.SW(true_value_locate(i)) info.SW(true_value_locate(i+1))],inter_x);
        info.SW(true_value_locate(i):true_value_locate(i+1)) = inter_approx;
    end 


    [~, post_pdf] = calpdf(samples, mu, sigma, ppi, ellipse, para_stru);
    [~, labels] = max(post_pdf, [], 2);
    %[~, labels] = max(pdfX, [], 2);
    samples_rec = samples;
    for k = 1:K
        temp = find(labels==k);
        len = length(temp);
        tempMU = repmat(mu(k,:),len,1); %mvtrnd(sigma_EMM(:,:,mm),1,len)+repmat(mu_EMM(mm,:),len,1);
        samples_rec(temp,:) = tempMU; 
    end
    final_mu = mu;
    final_sigma = sigma;
    final_pi = ppi;