% obtain optimal transport map for one-dimensional distributions
function [map, kpotential, sw] = optimal_map(pdf_p,pdf_e,range) % from dist p to dist e
    pdf_p = pdf_p + eps;
    pdf_e = pdf_e + eps;
    pdf_p = pdf_p./sum(pdf_p);
    pdf_e = pdf_e./sum(pdf_e);
    cdf_p = cumsum(pdf_p);
    cdf_e = cumsum(pdf_e);
    unix = 1:length(cdf_p);
    uniy = linspace(0,1,length(cdf_p));
    uni_cdf_p = interp1(cdf_p, unix, uniy, 'pchip');
    uni_cdf_e = interp1(cdf_e, unix, uniy, 'pchip');
    mapdiff = interp1(uni_cdf_p, uni_cdf_p-uni_cdf_e, unix, 'pchip');
    map = unix - mapdiff; 
    kpotential = cumsum(mapdiff)./length(mapdiff).*range;
    kpotential = kpotential - mean(kpotential);
    sw = sum((mapdiff./length(mapdiff).*range).^2.*pdf_p); % To rescale it into normal coordinates