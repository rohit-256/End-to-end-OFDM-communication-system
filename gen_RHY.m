function R_HYp = gen_RHY(h_var,N,index,p_indices) 
    n = 0 : length(h_var) - 1;

    temp = exp(1i*2*pi*(1/N)*(p_indices - index) .* transpose(n));
    R_HYp = temp .* transpose(h_var);
    R_HYp = sum(R_HYp,1);
end
