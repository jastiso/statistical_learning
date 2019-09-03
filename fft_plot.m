function [f,P1] = fft_plot(x,L,srate)
% makes an FFT plot
% input (x) should be a vector

Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
f = srate*(0:(L/2))/L;
plot(f,P1)

end

