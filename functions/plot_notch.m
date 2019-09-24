%% Plot notch filter
clear 

data = [];
srate = 1000;
f_data = filter_resp(data, 4000, srate);
img_dir = ['/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_img']

% filter out 60 Hz harmonics
[b,a] = butter(4, [59/(srate/2), 61/(srate/2)], 'stop');
data60 = filtfilt(b,a,f_data.trial{1}')';

[b,a] = butter(4, [119/(srate/2), 121/(srate/2)], 'stop');
data120 = filtfilt(b,a,f_data.trial{1}')';

[b,a] = butter(4, [179/(srate/2), 181/(srate/2)], 'stop');
data180 = filtfilt(b,a,f_data.trial{1}')';

figure(1); clf
subplot(2,2,1)
plot(data60(1,:)', 'linewidth', 2); hold on
plot(data120(1,:)', 'linewidth', 2);
plot(data180(1,:)', 'linewidth', 2);
legend([{'60 Hx'},{'120 Hz'},{'180 Hz'}])
title('Impulse Response')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,2)
plot(f_data.trial{1}(1,:)', 'k', 'linewidth', 2);
title('Impulse')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,3)
plot(data60(2,:)', 'linewidth', 2); hold on
plot(data120(2,:)', 'linewidth', 2);
plot(data180(2,:)', 'linewidth', 2);
legend([{'60 Hx'},{'120 Hz'},{'180 Hz'}])
title('Step Response')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,4)
plot(f_data.trial{1}(2,:)', 'k', 'linewidth', 2)
title('Step')
xlabel('Time')
ylabel('Amplitude')
saveas(gca, [img_dir, '/notch_filter_resp.png'], 'png')


%% Repeat for FFTs

L = 4000;
[freqs, fourier_impulse] = fft_plot(f_data.trial{1}(1,:), L, srate);
[~, fourier_step] = fft_plot(f_data.trial{1}(2,:), L, srate);
[~, fourier_imp60] = fft_plot(data60(1,:), L, srate);
[~, fourier_imp120] = fft_plot(data120(1,:), L, srate);
[~, fourier_imp180] = fft_plot(data180(1,:), L, srate);
[~, fourier_step60] = fft_plot(data60(2,:), L, srate);
[~, fourier_step120] = fft_plot(data120(2,:), L, srate);
[~, fourier_step180] = fft_plot(data180(2,:), L, srate);

figure(2); clf
subplot(1,2,1)
plot(freqs, fourier_imp60', 'linewidth', 2); hold on
plot(freqs, fourier_imp120', 'linewidth', 2);
plot(freqs, fourier_imp180', 'linewidth', 2);
plot(fourier_impulse', 'k', 'linewidth', 2);
legend([{'60 Hx'},{'120 Hz'},{'180 Hz'},{'No FIlter'}])
title('Impulse Response FFT')
xlabel('Freq (Hz)')
ylabel('Power')
subplot(1,2,2)
plot(fourier_step', 'k', 'linewidth', 2); hold on
plot(fourier_step60', 'linewidth', 2); 
plot(fourier_step120', 'linewidth', 2);
plot(fourier_step180', 'linewidth', 2);
legend([{'No filter'},{'60 Hx'},{'120 Hz'},{'180 Hz'}])
title('Step Response FFT')
xlabel('Freq (Hz)')
ylabel('Power')



