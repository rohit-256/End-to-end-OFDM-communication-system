function R_HH = gen_RHH(h_var,N,p_indices)
    
    R_HH = zeros(length(p_indices));
    for i = 1 : length(p_indices)
        R_HH(i,:) = gen_RHY(h_var,N,p_indices(i),p_indices) ;
    end
end