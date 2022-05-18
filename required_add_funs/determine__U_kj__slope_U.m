function [U_kj, slope_U] = determine__U_kj__slope_U(U_k,U_kp1,j,k,m,M_discr,N,h,U_discont)
    U_kj = U_k((j-1)*m+1:j*m);
    if j<M_discr
        U_kjp1 = U_k(j*m+1:(j+1)*m);
    elseif j==M_discr && k<N-1
        U_kjp1 = U_kp1(1:m);
    elseif j==M_discr && k==N-1
        if U_discont==1
            U_kjp1 = U_k((M_discr-1)*m+1:M_discr*m);
        else
            U_kjp1 = U_kp1(1:m);
        end
    end
    slope_U = (U_kjp1-U_kj)/h;
    slope_U2 = U_kjp1-U_kj;
    if U_discont==1
        slope_U  = 0*slope_U;
        slope_U2 = 0*slope_U2;
    end
end