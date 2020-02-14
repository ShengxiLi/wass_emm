% gradient of g(t)
function out_psi = psi_fun(in_t, in_ellip, in_para_stru, M)
%  sign_t = sign(in_t);
%  in_t = abs(in_t);
if strcmp(in_ellip, 'PVII')
    out_psi = -in_para_stru.s./(in_para_stru.v + in_t);
elseif strcmp(in_ellip, 'Kotz')
    out_psi = (in_para_stru.a - 1)./in_t - in_para_stru.b.*in_para_stru.s.*in_t.^(in_para_stru.s - 1);
elseif strcmp(in_ellip, 'Logi')
    out_psi = (exp(-in_t)-1)./(exp(-in_t)+1);
elseif strcmp(in_ellip, 'Hype')
    out_psi = (in_para_stru.lambda - (M/2))./(in_para_stru.a + in_t) - in_para_stru.v./2./((in_para_stru.a*in_para_stru.v + in_para_stru.v.*in_t).^0.5).*besselk((in_para_stru.lambda - (M/2) + 1),(in_para_stru.a*in_para_stru.v + in_para_stru.v.*in_t).^0.5)./besselk((in_para_stru.lambda - (M/2)),(in_para_stru.a*in_para_stru.v + in_para_stru.v.*in_t).^0.5);
elseif strcmp(in_ellip, 'PeII')
    out_psi = (in_para_stru.s - 1)./(in_t - 1);
end
% out_psi = out_psi.*sign_t;