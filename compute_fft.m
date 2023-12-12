function [f, fft_dB] = compute_fft(X, Fe)
    
    N = length(X);
    f = Fe/N*(0:(N/2));
    Y= fft(X);
    P2 = abs(Y/N);
    P1 = P2(1:fix(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    fft_dB = 20*log(P1);
end
