% index several types of EMMs
function [ellipse, parameters] = idx2dist(dist_idx, K, M)
% function to index an elliptical distribution
    ellipse = [];
    ellips_name = ['Kotz'; 'PVII'; 'Hype'; 'Logi'; 'PeII'];
    if (dist_idx == 1) % Gaussian distribution
        for k = 1:K
        ellipse = [ellipse; ellips_name(1,:)];
        parameters{k}.a = 1; parameters{k}.s = 1; parameters{k}.b = 0.5;
        end 
    elseif (dist_idx == 2) % Student-t (t=1)/Cauchy distribution
        for k = 1:K
        ellipse = [ellipse; ellips_name(2,:)];
        parameters{k}.v = 1; parameters{k}.s = (M+parameters{k}.v)/2;
        end     
    elseif (dist_idx == 3) % Logistic distribution
        for k = 1:K
            ellipse = [ellipse; ellips_name(4,:)];
            parameters{k}.non = 1.5; 
        end   
    elseif (dist_idx == 4) % Gamma distribution
        for k = 1:K
            ellipse = [ellipse; ellips_name(1,:)];
            parameters{k}.a = 2; parameters{k}.s = 1; parameters{k}.b = 0.5;
        end  
    end