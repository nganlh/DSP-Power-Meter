%% FFT THD Calculation with Full Spectrum Plot
clear; clc; close all;

%% --- Signal Settings ---
fs       = 2000;       % Sampling frequency (Hz)
f_signal = 50;         % Fundamental frequency (Hz)
duration = 0.2;        % Signal duration (s)
t        = (0:1/fs:duration-1/fs)';  % Column vector

% Voltage signal: fundamental + harmonics
V_amp = 220*sqrt(2);   % Peak amplitude for 220 Vrms
voltage = V_amp*sin(2*pi*f_signal*t) + ...
          10*sin(2*pi*3*f_signal*t) + ...
          5*sin(2*pi*5*f_signal*t);  % Add harmonics

N = length(t);          % number of samples

%% --- FFT Calculation ---
V_fft = fft(voltage)/N;          % FFT and normalize
V_mag = 2*abs(V_fft(1:N/2+1));   % Single-sided amplitude spectrum
f = (0:N/2)*fs/N;                % Frequency vector

% Find fundamental index
[~,fund_idx] = min(abs(f - f_signal));

% THD calculation (exclude fundamental)
V_thd = sqrt(sum(V_mag(2:end).^2) - V_mag(fund_idx)^2) / V_mag(fund_idx);
fprintf('Voltage THD (FFT) = %.3f %%\n', V_thd*100);

%% --- Plots ---

figure('Name','Voltage Signal and FFT','Units','normalized','Position',[0.05 0.05 0.9 0.7]);

% Time-domain waveform
subplot(3,1,1);
plot(t, voltage, 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Voltage [V]');
title('Voltage Signal');

% Zoomed harmonic magnitude (first 10 harmonics)
subplot(3,1,2);
stem(f(1:10), V_mag(1:10), 'filled','LineWidth',1.5); grid on;
xlabel('Frequency [Hz]'); ylabel('Amplitude [V]');
title('Harmonic Magnitudes (FFT)');
xlim([0 10*f_signal]);

% Full FFT spectrum
subplot(3,1,3);
plot(f, V_mag, 'LineWidth', 1.5); grid on;
xlabel('Frequency [Hz]'); ylabel('Amplitude [V]');
title('Full FFT Spectrum');
xlim([0 fs/2]);  % Up to Nyquist frequency
