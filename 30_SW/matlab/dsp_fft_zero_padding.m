clc; clear; close all;

% Parameters
fs = 7200;        % Sampling frequency (Hz)
f_sig = 60;      % Signal frequency (Hz)
N = 720;         % Number of real samples (not power of 2)

t = (0:N-1)/fs;  % Time vector

% Generate 50 Hz sine wave
x = sin(2*pi*f_sig*t);

%% FFT with 720 samples
X1 = fft(x, N);              % FFT without zero padding
f1 = (0:N-1)*(fs/N);         % Frequency axis
mag1 = abs(X1)/N;            % Magnitude spectrum (normalize)

%% FFT with zero padding to 1024
N_pad = 8192;
X2 = fft(x, N_pad);          % FFT with zero padding
f2 = (0:N_pad-1)*(fs/N_pad); % Frequency axis
mag2 = abs(X2)/N;            % Note: normalize by original N

%% Plot results
figure;

subplot(2,1,1);
stem(f1, mag1, 'b', 'Marker', 'none'); hold on;
xlim([0 200]); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT with 720 samples (no padding)');

subplot(2,1,2);
plot(f2, mag2, 'r'); hold on;
xlim([0 200]); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT with zero padding to 1024');

