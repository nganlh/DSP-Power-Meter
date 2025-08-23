%% Goertzel Example for THD Calculation
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

%% --- Goertzel Setup ---
num_harm = 10;                               % fundamental + 9 harmonics
k = round((f_signal*(1:num_harm))*N/fs);     % bin indices for harmonics

V_goer = zeros(num_harm,1);

% Compute Goertzel for each harmonic
for h = 1:num_harm
    V_goer(h) = abs(goertzel(voltage, k(h)));
end

% Normalize amplitude (single-sided)
V_mag = (2/N) * V_goer;

%% --- THD Calculation ---
V_thd = sqrt(sum(V_mag(2:end).^2)) / V_mag(1);
fprintf('Voltage THD = %.3f %%\n', V_thd*100);

%% --- Plots ---
figure('Name','Voltage Signal and Harmonics','Units','normalized','Position',[0.05 0.05 0.9 0.6]);

% Time-domain waveform
subplot(2,1,1);
plot(t, voltage, 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Voltage [V]');
title('Voltage Signal');

% Harmonic magnitude plot
subplot(2,1,2);
stem(1:num_harm, V_mag, 'filled','LineWidth',1.5); grid on;
xlabel('Harmonic Order'); ylabel('Amplitude [V]');
title('Harmonic Magnitudes (Goertzel)');
xlim([0 num_harm+1]);
