% ==========================================
% AC Power Analysis using FFT (RMS, Power)
% ==========================================
clear; clc; close all;

% --- Parameters ---
f_signal = 50;          % Fundamental frequency (Hz)
fs = 7200;              % Sampling frequency (Hz)
N_signal = 720;         % Actual number of samples captured
N_fft = 1024;           % FFT size (zero padded)
T = N_signal / fs;      % Duration (s)

V_amp = 220 * sqrt(2);  % Voltage amplitude (for Vrms = 220 V)
I_amp = 100 * sqrt(2);  % Current amplitude (for Irms = 100 A)
phi = deg2rad(30);      % Phase shift (current lags voltage)

% --- Time vector (only for captured samples) ---
t = (0:N_signal-1)/fs;

% --- Generate fundamental signals ---
voltage = V_amp * sin(2*pi*f_signal*t);
current = I_amp * sin(2*pi*f_signal*t - phi);

% --- Add harmonics (optional) ---
% V3 = 0.2 * V_amp;   % 20% of fundamental
% V5 = 0.1 * V_amp;   % 10% of fundamental
% I3 = 0.15 * I_amp;  % 15% of fundamental
% I5 = 0.05 * I_amp;  % 5% of fundamental
% voltage = voltage + V3*sin(2*pi*(3*f_signal)*t) + V5*sin(2*pi*(5*f_signal)*t);
% current = current + I3*sin(2*pi*(3*f_signal)*t - phi) + I5*sin(2*pi*(5*f_signal)*t - phi);

% --- Zero padding to 1024 samples ---
voltage_padded = [voltage, zeros(1, N_fft - N_signal)];
current_padded = [current, zeros(1, N_fft - N_signal)];

% --- FFT calculation ---
V_fft_c = fft(voltage_padded)/N_signal;  % normalize by actual N_signal
I_fft_c = fft(current_padded)/N_signal;
f = (0:N_fft/2)*fs/N_fft;               % frequency vector (up to Nyquist)

% --- FFT magnitude & phase ---
V_fft = 2*abs(V_fft_c(1:N_fft/2+1));
I_fft = 2*abs(I_fft_c(1:N_fft/2+1));
V_phase = angle(V_fft_c(1:N_fft/2+1));
I_phase = angle(I_fft_c(1:N_fft/2+1));

% --- Compute RMS from FFT ---
% RMS = sqrt(sum((magnitude/sqrt(2)).^2))
Vrms_fft = sqrt(sum((V_fft(2:end)/sqrt(2)).^2));  % skip DC
Irms_fft = sqrt(sum((I_fft(2:end)/sqrt(2)).^2));  % skip DC

% --- Compute Real Power (from all harmonics) ---
P_fft = 0;
for k = 2:length(V_fft)
    P_fft = P_fft + (V_fft(k)/2) * (I_fft(k)/2) * cos(V_phase(k) - I_phase(k));
end

% --- Apparent Power, Power Factor, Reactive Power ---
S_fft = Vrms_fft * Irms_fft;
PF_fft = P_fft / S_fft;
Q_fft = sqrt(max(S_fft^2 - P_fft^2, 0));

% --- Print results ---
fprintf("==== Frequency-Domain Power Results (Zero Padding) ====\n");
fprintf("Samples (signal): %d, FFT size: %d\n", N_signal, N_fft);
fprintf("Vrms (FFT): %.3f V\n", Vrms_fft);
fprintf("Irms (FFT): %.3f A\n", Irms_fft);
fprintf("Real Power (FFT): %.3f W\n", P_fft);
fprintf("Apparent Power (FFT): %.3f VA\n", S_fft);
fprintf("Reactive Power (FFT): %.3f var\n", Q_fft);
fprintf("Power Factor (FFT): %.3f\n", PF_fft);

% --- Plot FFT Magnitude ---
figure;
subplot(2,1,1);
stem(f, V_fft, 'b', 'LineWidth', 1.2); hold on;
stem(f, I_fft, 'r', 'LineWidth', 1.2);
xlim([0 500]);
legend('Voltage FFT','Current FFT');
title('FFT Magnitude Spectrum (Zero-Padded)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
grid on;

% --- Plot phase difference per frequency ---
subplot(2,1,2);
plot(f, rad2deg(V_phase - I_phase), 'k', 'LineWidth', 1.2);
xlim([0 500]);
title('Phase Difference (Voltage - Current)');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;

% --- Frequency Resolution Comparison ---
df_no_pad = fs / N_signal;
df_pad = fs / N_fft;
fprintf("\nFrequency resolution without zero padding: %.2f Hz\n", df_no_pad);
fprintf("Frequency resolution with zero padding:    %.2f Hz\n", df_pad);
