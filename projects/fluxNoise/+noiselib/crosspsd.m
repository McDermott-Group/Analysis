function [psd] = crosspsd(z1,z2,Fs)
%Cross PSD

NumP = length(z1);
fft_seq1 = fft(z1);
fft_seq2 = conj(fft(z2));
psd = (1/(Fs*NumP))*(fft_seq1.*fft_seq2);
psd(2:(NumP/2)) = 2*psd(2:(NumP/2));

end