% ==========================================
% FFT-based Vrms/Irms/Power using zero padding
% Correct normalization using Parseval
% ==========================================
clear; clc; close all;

% Parameters
f_signal  = 50;          % Hz
fs        = 7200;        % Hz
N_signal  = 720;         % actual captured samples
N_fft     = 1024;        % FFT length (zero padded)
t         = (0:N_signal-1)/fs;

V_amp = 220 * sqrt(2);   % peak for 220 Vrms
I_amp = 100 * sqrt(2);   % peak for 100 Arms
phi   = deg2rad(30);     % current lags

% Generate signals (720 samples)
voltage = V_amp * sin(2*pi*f_signal*t);
current = I_amp * sin(2*pi*f_signal*t - phi);

% % (optional) add harmonics if you want:
% V3 = 0.2 * V_amp; V5 = 0.1 * V_amp;
% I3 = 0.15 * I_amp; I5 = 0.05 * I_amp;
% voltage = voltage + V3*sin(2*pi*(3*f_signal)*t) + V5*sin(2*pi*(5*f_signal)*t);
% current = current + I3*sin(2*pi*(3*f_signal)*t - phi) + I5*sin(2*pi*(5*f_signal)*t - phi);

% Zero-pad to N_fft
voltage_p = [voltage, zeros(1, N_fft - N_signal)];
current_p = [current, zeros(1, N_fft - N_signal)];

% Compute complex FFT (no /N normalization here)
Vc = fft(voltage_p, N_fft);
Ic = fft(current_p, N_fft);

% Frequency axis (up to Nyquist for plotting)
f = (0:(N_fft/2)) * (fs / N_fft);

% Single-sided magnitudes for plotting (optional)
V_mag_ss = 2*abs(Vc(1:N_fft/2+1))/N_signal;   % divide by N_signal to get amplitudes
I_mag_ss = 2*abs(Ic(1:N_fft/2+1))/N_signal;

% ----- CORRECT Vrms, Irms using Parseval -----
% Vrms^2 = (1/(N_signal*N_fft)) * sum_k |Vc(k)|^2
Vrms_fft = sqrt( sum(abs(Vc).^2) / (N_signal * N_fft) );
Irms_fft = sqrt( sum(abs(Ic).^2) / (N_signal * N_fft) );

% ----- CORRECT real power from FFT -----
% P = (1/(N_signal*N_fft)) * Re{ sum_k Vc(k) * conj(Ic(k)) }
P_fft = real( sum( Vc .* conj(Ic) ) ) / (N_signal * N_fft);

% Apparent and reactive
S_fft = Vrms_fft * Irms_fft;
Q_fft = sqrt(max(S_fft^2 - P_fft^2, 0));
PF_fft = P_fft / S_fft;

% Print results (compare to time-domain)
Vrms_time = sqrt(mean(voltage.^2));
Irms_time = sqrt(mean(current.^2));
P_time = mean(voltage .* current);
S_time = Vrms_time * Irms_time;
PF_time = P_time / S_time;

fprintf('Time-domain: Vrms=%.6f V, Irms=%.6f A, P=%.6f W, PF=%.6f\n', ...
    Vrms_time, Irms_time, P_time, PF_time);
fprintf('FFT-domain : Vrms=%.6f V, Irms=%.6f A, P=%.6f W, PF=%.6f\n', ...
    Vrms_fft, Irms_fft, P_fft, PF_fft);

% ----- Optional: THD from single-sided magnitudes (use proper scaling) -----
% compute single-sided amplitudes in peak units (already computed V_mag_ss)
% fundamental bin index (nearest)
[~, fund_bin] = min(abs(f - f_signal));
V1 = V_mag_ss(fund_bin);
I1 = I_mag_ss(fund_bin);
V_thd = sqrt(sum(V_mag_ss(fund_bin*2:end).^2)) / V1; % approximate
I_thd = sqrt(sum(I_mag_ss(fund_bin*2:end).^2)) / I1;

fprintf('Voltage THD (approx) = %.4f %%\n', V_thd*100);
fprintf('Current THD (approx) = %.4f %%\n', I_thd*100);

% ----- Plots -----
figure('Name','FFT mag (single-sided)');
subplot(2,1,1);
plot(f, V_mag_ss, 'b'); xlim([0 500]); grid on;
title('Voltage single-sided amplitude (zero-padded)');
xlabel('Hz'); ylabel('Amplitude (peak)');

subplot(2,1,2);
plot(f, I_mag_ss, 'r'); xlim([0 500]); grid on;
title('Current single-sided amplitude (zero-padded)');
xlabel('Hz'); ylabel('Amplitude (peak)');
