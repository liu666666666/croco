function [s_data_v] = fftsmooth(data_v,freq_n)
% use fft low pass filter to smoothen out a signal
%data_v=data_v-mean(data_v);
%data_v=data_v.*tukeywin(length(data_v),1)';
fft_data_v = fft(data_v);
s_fft_data_v = zeros(1,length(data_v));
s_fft_data_v(1:freq_n) = fft_data_v(1:freq_n);
s_fft_data_v(end-freq_n:end) = fft_data_v(end-freq_n:end);
s_data_v = real(ifft(s_fft_data_v));

