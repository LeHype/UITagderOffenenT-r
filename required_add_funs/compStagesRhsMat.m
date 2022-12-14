function stages_rhs_mat = compStagesRhsMat(IntMethod,stages_mat,rhs,t,x,u,u_slope,h,butcher,rhsParams)

A = butcher.A;
c = butcher.c;
s = length(c);

stages_rhs_mat = 0*casadi.SX.sym('stages_rhs_mat',rhs.size1_out(0),s);

for ii=1:s
    stage_sum_i = 0*x;
    switch IntMethod
        case 'explicit'
            for jj=1:ii-1
                stage_sum_i = stage_sum_i + A(ii,jj)*stages_rhs_mat(:,jj);
            end  
        case 'implicit'
            for jj=1:s
                stage_sum_i = stage_sum_i + A(ii,jj)*stages_mat(:,jj);
            end    
    end
    if ~isempty(rhsParams)
        stages_rhs_mat(:,ii) = rhs(t + c(ii)*h,...
                                   x + h*stage_sum_i,...
                                   u + c(ii)*h*u_slope,...
                                   rhsParams{:});
    else
        stages_rhs_mat(:,ii) = rhs(t + c(ii)*h,...
                                   x + h*stage_sum_i,...
                                   u + c(ii)*h*u_slope);
    end
end


% k_mat = 0*repmat(Y,1,s);    
% for iii=1:s
%     stages_sum_i = 0*Y;
%     for jjj=1:iii-1
%         stages_sum_i = stages_sum_i + butcher.A(iii,jjj)*k_mat(:,jjj);
%     end
%     k_mat(:,iii) = rhs_Mayer(t+butcher.c(iii)*h, Y + h*stages_sum_i, U_0j+butcher.c(iii)*h*slope_U,w_L);
% end

end