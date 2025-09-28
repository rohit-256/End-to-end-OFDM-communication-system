function constellation = MPSK(M)
    
    k = 1 : (M);
    constellation = exp(2i*pi*(k-1)/M);
    
     
end