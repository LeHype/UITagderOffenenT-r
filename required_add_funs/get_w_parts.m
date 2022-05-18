function w_parts = get_w_parts(w,n,m,M,use_FOH)
    if use_FOH==1
        length_w_extra = m;
    else
        length_w_extra = 0;
    end
    length_w_time = 1;
    w_middle_SE_inds = [(M*m+n)+1,...
                        length(w)-length_w_extra-(m+n)-length_w_time];
    w_parts.start  = w(1 : M*m+n);
    w_parts.middle = w(w_middle_SE_inds(1) : w_middle_SE_inds(2));
    if use_FOH==1
        w_parts.extra = w(w_middle_SE_inds(2)+1 : w_middle_SE_inds(2)+length_w_extra);
    else
        w_parts.extra = [];
    end
    w_parts.end    = w(end-(m+n)+1-length_w_time : end-length_w_time);
    w_parts.time   = w(end-length_w_time+1:end);
end