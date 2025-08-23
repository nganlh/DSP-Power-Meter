% ==========================================
% AC Power Analysis in MATLAB + THD plots
% ==========================================

clear; clc; close all;

% --- Parameters ---
f_signal = 50;        % Fundamental frequency (Hz)
fs = 5000;            % Sampling frequency (Hz)
T = 0.2;              % Simulation time (s)

V_amp = 220 * sqrt(2);    % Voltage amplitude (Volts)
I_amp = 100 * sqrt(2);      % Current amplitude (Amps)
phi = deg2rad(30);        % Phase shift of current (radians, e.g., 30° lag)

% --- Time vector ---
t = 0:1/fs:T;

% --- Generate signals ---
voltage = V_amp * sin(2*pi*f_signal*t);
current = I_amp * sin(2*pi*f_signal*t - phi);

% Add harmonics (example: 3rd and 5th)
V3 = 0.2 * V_amp;   % 20% amplitude of fundamental
V5 = 0.1 * V_amp;   % 10% amplitude
I3 = 0.15 * I_amp;  % 15% amplitude
I5 = 0.05 * I_amp;  % 5% amplitude
voltage = voltage + V3*sin(2*pi*(3*f_signal)*t) + V5*sin(2*pi*(5*f_signal)*t);
current = current + I3*sin(2*pi*(3*f_signal)*t - phi) + I5*sin(2*pi*(5*f_signal)*t - phi);

% --- RMS calculations ---
Vrms = sqrt(mean(voltage.^2));
Irms = sqrt(mean(current.^2));

% --- Power calculations ---
P = mean(voltage .* current);   % Real power (W)
S = Vrms * Irms;                % Apparent power (VA)
Q = sqrt(max(S^2 - P^2, 0));    % Reactive power (VAR)
PF = P / S;                     % Power factor

% --- THD calculations (using FFT) ---
N = length(t);
V_fft = abs(fft(voltage))/N;
I_fft = abs(fft(current))/N;

V_fft = V_fft(1:N/2+1);   % Keep positive frequencies
I_fft = I_fft(1:N/2+1);

f = (0:N/2)*fs/N;         % Frequency vector

% Find fundamental index
[~,fund_idx] = min(abs(f - f_signal));

% THD = sqrt(sum(harmonics^2)) / fundamental
V_thd = sqrt(sum(V_fft(fund_idx*2:end).^2)) / V_fft(fund_idx);
I_thd = sqrt(sum(I_fft(fund_idx*2:end).^2)) / I_fft(fund_idx);

% --- Results ---
fprintf("==== Power Analysis Results ====\n");
fprintf("N: %d\n", N);
fprintf("Vrms: %.2f V\n", Vrms);
fprintf("Irms: %.2f A\n", Irms);
fprintf("Apparent Power (S): %.2f VA\n", S);
fprintf("Real Power (P): %.2f W\n", P);
fprintf("Reactive Power (Q): %.2f var\n", Q);
fprintf("Power Factor (PF): %.3f\n", PF);
fprintf("Voltage THD: %.3f %%\n", V_thd*100);
fprintf("Current THD: %.3f %%\n", I_thd*100);

% --- Plot signals ---
figure;
subplot(2,2,1);
plot(t, voltage, 'b', 'LineWidth', 1.2); hold on;
plot(t, current, 'r', 'LineWidth', 1.2);
legend('Voltage','Current');
title('AC Voltage & Current');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(2,2,2);
stem(f, V_fft, 'b'); hold on;
stem(f, I_fft, 'r');
xlim([0 500]);
legend('Voltage Spectrum','Current Spectrum');
title('FFT Spectrum');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
grid on;

% --- Voltage & Current Harmonic Spectrum (Separate Subplots) ---
harmonics = 1:10; % show first 10 harmonics
V_harm = zeros(size(harmonics));
I_harm = zeros(size(harmonics));

for k = harmonics
    idx = fund_idx * k;
    if idx <= length(V_fft)
        V_harm(k) = V_fft(idx);
    end
    if idx <= length(I_fft)
        I_harm(k) = I_fft(idx);
    end
end

% Voltage subplot
subplot(2,2,3);
bar(harmonics, V_harm, 'b');
title(sprintf('Voltage Harmonic Spectrum (THD = %.2f%%)', V_thd*100));
xlabel('Harmonic Number');
ylabel('Amplitude');
grid on;

% Current subplot
subplot(2,2,4);
bar(harmonics, I_harm, 'r');
title(sprintf('Current Harmonic Spectrum (THD = %.2f%%)', I_thd*100));
xlabel('Harmonic Number');
ylabel('Amplitude');
grid on;
