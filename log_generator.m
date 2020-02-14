% output c_m and g(t) for each type of EMM, in log format to facilitate
% small value calculations
function out_log_generator = log_generator(in_t, in_ellip, in_para_stru, M)
%  in_t = abs(in_t);
if strcmp(in_ellip, 'PVII')
    CM = log((pi*in_para_stru.v)^(-M/2)*gamma(in_para_stru.s)/gamma(in_para_stru.s - (M/2)));
    out_log_generator = CM + (-in_para_stru.s).*log(1+ in_t/in_para_stru.v);
elseif strcmp(in_ellip, 'Kotz')
    CM = log(gamma(M/2)*in_para_stru.s*in_para_stru.b^((2*in_para_stru.a+M-2)/2/in_para_stru.s)/gamma((2*in_para_stru.a+M-2)/2/in_para_stru.s)/pi^(M/2));
    out_log_generator = CM + (in_para_stru.a - 1).*log(in_t) - in_para_stru.b.*in_t.^in_para_stru.s;
elseif strcmp(in_ellip, 'Logi')
    out_log_generator = -in_t - 2.*log(1+exp(-in_t));
elseif strcmp(in_ellip, 'Hype')
    CM = log((in_para_stru.v/in_para_stru.a)^(in_para_stru.lambda/2)/(2*pi)^(M/2)/besselk(in_para_stru.lambda,(in_para_stru.v*in_para_stru.a)^0.5));
    out_log_generator = CM + (in_para_stru.lambda/2 - M/4).*log(in_para_stru.a/in_para_stru.v + in_t./in_para_stru.v) + log(besselk((in_para_stru.lambda - M/2),(in_para_stru.a*in_para_stru.v + in_t.*in_para_stru.v).^0.5));
elseif strcmp(in_ellip, 'PeII')
    CM = log(gamma(M/2+in_para_stru.s)/pi^(M/2)/gamma(in_para_stru.s));
    out_log_generator = CM + (in_para_stru.s - 1).*log(1 - in_t);
end

if sum(isnan(out_log_generator))
    error('log_generator is nan');
end

if sum(isinf(out_log_generator))
    error('log_generator is inf');
end