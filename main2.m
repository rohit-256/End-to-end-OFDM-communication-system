
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%indices start from 0 (0 to N-1)
%pilot indices are 0, S-1, 2s-1,.....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = main2(N,S,CP,L,index,M)
    tic
%     N = str2double(N)
%     S = str2double(S)
%     CP = str2double(CP)
%     L = str2double(L)
%     index = str2double(index)
%     M = str2double(M)
    
    %pilot locations
    index = index + 1;
    Np = N/S;
    p_indices = [1,(1 : Np-1)*S];

    %method2 and 3 indices
    if ismember(index,p_indices)
        p_indices2 = index;
        i1 = find(p_indices == index);
    elseif index > (Np-1)*S
        p_indices2 = (Np-1)*S;
        i1 = Np;
        
    else
        i1 = find(p_indices<index, 1, 'last');
        i2 = i1 + 1;
        p_indices2 = p_indices([i1,i2]);
    end

    
    sample_size = 1e5;
    snr_db = 0 : 2 : 40;
    snr = 10.^(snr_db/10);
    
    



    constellation = MPSK(M);
    
    tx_vector = constellation(randi(M,1,sample_size*N));
    tx_vector = transpose(reshape(tx_vector,N,sample_size));
    tx_vector(:,p_indices) = 1;

    %ifft
    tx_idft = sqrt(N)*ifft(tx_vector,N,2);
    
    %add cp

    tx_cp = [tx_idft(:,N-CP+1:N) tx_idft];


    %%%%Rayleigh Channel%%%%
    m = 1:L;
    h_var = ((0.5).^m)/(1-(0.5)^L);

    h = sqrt(1/2)*randn(sample_size,L) + sqrt(1/2)*1i*randn(sample_size,L);
    h = h .* sqrt(h_var);
    H = fft(h,N,2);
    
    %MMSE matrices
    R_HYp = gen_RHY(h_var,N,index,p_indices);
    R_HH = gen_RHH(h_var,N,p_indices);

    %method3 R_HY
    if index > (Np-1)*S
        R_HYp3 = R_HH(i1,i1);
    
    else
        R_HYp3 = R_HH(i1,i1);
        R_HYp4 = R_HH(i2,i2);

    end

    %transmission and decoding


    Pe = zeros(1,length(snr));
    Pe1 = zeros(1,length(snr));
    Pe2 = zeros(1,length(snr));
    Pe3 = zeros(1,length(snr));


    MSE1 = zeros(1,length(snr));
    MSE2 = zeros(1,length(snr));
    MSE3 = zeros(1,length(snr));
    n = sqrt(1/2)*randn(sample_size,N + CP + L-1) + sqrt(1/2)*1i*randn(sample_size,N + CP + L-1);

    for i = 1: length(snr)
        n_i = sqrt(1/snr(i)) * n;
        
        rx = block_conv(h,tx_cp) + n_i;
        
        %Remove CP
        rx = rx(:,CP+1 : CP + N);

        %dft
        rx_dft = sqrt(1/N) * fft(rx,N,2);
        

        %Estimation
        %Method 1
        rx_pilots = rx_dft(:,p_indices);
        R_YpYp = R_HH + eye(Np)*(1/snr(i));
        
        Hi_hat = transpose(R_HYp /(R_YpYp) * transpose(rx_pilots));
        
        MSE1(i) = mean((abs(H(:,index) - Hi_hat)).^2);
        
        %method2
        if (numel(p_indices2) == 2)
            rx_pilots2 = rx_dft(:,p_indices2);
            R_YpYp2 = R_YpYp([i1,i2],[i1,i2]);
            R_HYp2 = R_HYp([i1,i2]);
            
            Hi_hat2 = transpose(R_HYp2 /(R_YpYp2) * transpose(rx_pilots2));
        else
            rx_pilots2 = rx_dft(:,p_indices2);
            R_YpYp2 = R_YpYp(i1,i1);
            R_HYp2 = R_HYp(i1);
            
%             size(R_HYp2)
%             size(R_YpYp2)
%             size(rx_pilots2)
            Hi_hat2 = transpose(R_HYp2 /(R_YpYp2) * transpose(rx_pilots2));
        end
        MSE2(i) = mean((abs(H(:,index) - Hi_hat2)).^2);


        %method3
        if (numel(p_indices2) == 2)
            rx_pilots_first = rx_dft(:,p_indices2(1));
            R_YpYp3 = R_YpYp2(1,1);
            
            Hi_hat_first = transpose(R_HYp3 /(R_YpYp3) * transpose(rx_pilots_first));
    
            rx_pilots_second = rx_dft(:,p_indices2(2));
            R_YpYp4 = R_YpYp2(2,2);
            Hi_hat_second = transpose(R_HYp4  /(R_YpYp4) * transpose(rx_pilots_second));
            
            m = (Hi_hat_second - Hi_hat_first)/(p_indices2(2) - p_indices2(1));
            Hi_hat3 = m*(index - p_indices2(1)) + Hi_hat_first;
        else
            rx_pilots_first = rx_dft(:,(Np - 1)*S);
            R_YpYp3 = R_YpYp(Np,Np);
            
            Hi_hat3 = transpose(R_HYp3 /(R_YpYp3) * transpose(rx_pilots_first));
        end
        MSE3(i) = mean((abs(H(:,index) - Hi_hat3)).^2);



        %normalization and SER
        %Perfect CSI
        rx_norm = rx_dft./H;
        
        rx_dec = psk_dec(rx_norm,constellation,M);
        
        
        Pe(i) = (nnz(rx_dec(:,index) ~= tx_vector(:,index)))/sample_size;
        
        %method 1
        rx_norm = rx_dft./Hi_hat;
        
        rx_dec = psk_dec(rx_norm,constellation,M);
        
        
        Pe1(i) = (nnz(rx_dec(:,index) ~= tx_vector(:,index)))/sample_size;
        
        %method 2
        rx_norm = rx_dft./Hi_hat2;
        
        rx_dec = psk_dec(rx_norm,constellation,M);
        
        
        Pe2(i) = (nnz(rx_dec(:,index) ~= tx_vector(:,index)))/sample_size;

        %method 3
        rx_norm = rx_dft./Hi_hat3;
        
        rx_dec = psk_dec(rx_norm,constellation,M);
        
        
        Pe3(i) = (nnz(rx_dec(:,index) ~= tx_vector(:,index)))/sample_size;

    end
figure;
semilogy(snr_db,Pe)
hold on
semilogy(snr_db,Pe1)
hold on
semilogy(snr_db,Pe2)
hold on
semilogy(snr_db,Pe3)
hold off
legend("Perfect CSI","all pilots","nearest 2 neighbours","linear interpolation")
title(sprintf("SER for subcarrier index %d",index-1));
xlabel("SNR")
ylabel("SER")
grid on

figure;
semilogy(snr_db,MSE1)
hold on
semilogy(snr_db,MSE2)
hold on
semilogy(snr_db,MSE3)
hold off
grid on
legend("all pilots","nearest 2 neighbours","linear interpolation")
title(sprintf("MSE for subcarrier index %d",index-1));
xlabel("SNR")
ylabel("SER")

toc
end