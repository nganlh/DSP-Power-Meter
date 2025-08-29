clear; clc; close all;

% --- Parameters ---
fs = 10000;             % Sampling frequency (Hz)
f_signal = 50;          % Fundamental frequency (Hz)
T = 10;                  % Duration (s)
t = 0:1/fs:T-1/fs;      % Time vector

V_amp = 220*sqrt(2);    % Voltage amplitude
I_amp = 10*sqrt(2);     % Current amplitude
phi   = 30*pi/180;      % Phase shift (30 degrees)

% --- Generate signals with harmonics ---
% Fundamental
voltage = V_amp * sin(2*pi*f_signal*t);
current = I_amp * sin(2*pi*f_signal*t - phi);

% Add 3rd & 5th harmonics to current (distortion)
% I3 = 0.15*I_amp; 
% I5 = 0.10*I_amp;
% current = current + I3*sin(2*pi*3*f_signal*t - phi) + ...
%                   I5*sin(2*pi*5*f_signal*t - phi);

% --- FFT Analysis ---
N = length(t);
f = (0:N-1)*(fs/N);

% Voltage FFT
V_fft = fft(voltage)/N;
I_fft = fft(current)/N;

% Magnitude
V_mag = abs(V_fft(1:N/2));
I_mag = abs(I_fft(1:N/2));
freqs = f(1:N/2);

% Find fundamental index
[~,fund_idx] = min(abs(freqs - f_signal));

% Extract fundamental phasors
V1 = V_fft(fund_idx+1);   % FFT is 0-based in MATLAB indexing
I1 = I_fft(fund_idx+1);

% Phase shift between fundamental voltage & current
phi_meas = angle(V1) - angle(I1);

% --- THD calculation for current ---
% RMS of harmonics
harmonics = [3 5]; % could extend
harm_idx = round(harmonics*f_signal*N/fs) + 1;
I_harm_rms = sqrt(sum(abs(I_fft(harm_idx)).^2))*sqrt(2); % convert to RMS

% Fundamental RMS
I1_rms = abs(I1)*sqrt(2);

THD_I = I_harm_rms / I1_rms;

% --- Power Factor ---
PF = cos(phi_meas) / sqrt(1 + THD_I^2);

% --- Results ---
fprintf('Measured phase shift = %.2f deg\n', phi_meas*180/pi);
fprintf('THD(I) = %.2f %%\n', THD_I*100);
fprintf('Power Factor = %.4f\n', PF);

% --- Plot ---
figure;
subplot(2,1,1);
plot(t(1:1000), voltage(1:1000), 'b', t(1:1000), current(1:1000), 'r');
legend('Voltage','Current');
xlabel('Time (s)'); ylabel('Amplitude');
title('Waveforms');

subplot(2,1,2);
stem(freqs, I_mag);
xlim([0 500]);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Current Spectrum');
