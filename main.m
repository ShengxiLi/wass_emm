clear all
clc
dbstop if error
warning off

% load data
load('flower_shaped_data.mat')

% Specify parameters
M = 2; K = 8; N = size(samples,1); 
% Specify the Gaussian distribution
[ellipse, parameters] = idx2dist(1,K,M);


% Initialisation, parameter formats: mu: KxM; sigma: MxMxK
[ini_mu, ini_sigma, ini_pi] = initial_data(samples, K, 0); 


% Optimisation cfg
cfg.maxite = 2000; % maximal iterations
cfg.cp = 50; % check status per 50 epochs
cfg.lr = 0.01; % learning rate

% Optimise parameters
[mu, sigma, ppi, info, labels, ~] = mixture_solver(samples, ini_mu, ini_sigma, ini_pi, ellipse, parameters, cfg);
figure(1)
for k = 1:K
    temp = find(labels==k);
    part_samples = samples(temp,:);
    plot(part_samples(:,1), part_samples(:,2), '.'); hold on
end
for k = 1:K
    error_ellipse_adv(sigma(1:2,1:2,k), mu(k,1:2), 0.95, 'style', 'k'); hold on
end
xlim([-25,25]); ylim([-25,25]);
set(gca, 'fontsize', 30, 'fontname', 'Times New Roman', 'fontweight', 'Bold');

LL_ALL{1} = info.LL; % you can also use LL_ALL{2} = info2.LL to plot multiple algorithms
name_file{1} = 'Our Method'; % you can also use name_file{2} = 'Compared Methods' to plot multiple algorithms
% the same usage as above, in which you can easily extend to multiple plots
% by adding cells
MU_ALL{1} = info.mu; SIG_ALL{1} = info.sigma; THETA_ALL{1} = info.theta; SW_ALL{1} = info.SW; 
dynamic_pics(samples, name_file, MU_ALL, SIG_ALL, THETA_ALL, SW_ALL, LL_ALL)