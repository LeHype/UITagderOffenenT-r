function dpp = diff_pp(pp)
    if pp.order>1
        dpp_coeffs = ((size(pp.coefs,2)-1):-1:0) .* pp.coefs;
        dpp_coeffs = dpp_coeffs(:,1:end-1);
        dpp        = pp;
        dpp.coefs  = dpp_coeffs;
        dpp.order  = pp.order-1;
    else
        error('Order of polynomial is 1.')
    end
end