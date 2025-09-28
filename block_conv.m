function [out] = block_conv(h,x)
    %x = sample size x (N + l-1)
    %h = sample size x l
    sample_size = size(h,1);
    n = size(h,2);
    m = size(x,2);
    h = [h,zeros(sample_size,m-1)];
    x = [x,zeros(sample_size,n-1)];
    H = fft(h,n+m-1,2);
    X = fft(x,n+m-1,2);
    OUT = H.*X;
    out = ifft(OUT,n+m-1,2);

end