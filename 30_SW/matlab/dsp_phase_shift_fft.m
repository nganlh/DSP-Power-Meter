clear; clc; close all;

%% --- Parameters ---
fs = 7200;             % Sampling frequency (Hz)
f_signal = 60;         % Fundamental frequency (Hz)
N = 720;              % FFT size (power of 2)
t = (0:N-1)/fs;        % Time vector

%% --- Generate signals ---
V_amp = 220 * sqrt(2);     % Peak voltage (for 220 Vrms)
I_amp = 100 * sqrt(2);     % Peak current (for 100 Arms)
phi = deg2rad(-30.1234);         % Current lags voltage by 30 degrees

v = V_amp * sin(2*pi*f_signal*t);
%i = I_amp * sin(2*pi*f_signal*t - phi);
i_pure = I_amp * sin(2*pi*f_signal*t - phi);

% Generate PWM waveform
f_pwm = 20e3;                 % PWM frequency
pwm_carrier = square(2*pi*f_pwm*t, 50);   % 50% duty PWM: +1 / -1
pwm_carrier = (pwm_carrier + 1) / 2;  % Convert from [-1,1] → [0,1]

% Apply PWM to current
i = i_pure .* pwm_carrier;   % Modulated current


%% --- Remove DC offset (important for real-world data) ---
v = v - mean(v);
i = i - mean(i);

%% --- Apply FFT ---
V_fft = fft(v) / N;
I_fft = fft(i) / N;

%% --- Compute single-sided spectra ---
V_mag = 2 * abs(V_fft(1:N/2+1));
I_mag = 2 * abs(I_fft(1:N/2+1));
V_phase = angle(V_fft(1:N/2+1));
I_phase = angle(I_fft(1:N/2+1));
f = (0:N/2) * fs / N;

%% --- Find fundamental bin (50 Hz) ---
[~, k_fund] = min(abs(f - f_signal));

%% --- Calculate phase shift ---
phase_v = V_phase(k_fund);
phase_i = I_phase(k_fund);
phase_diff = phase_v - phase_i;       % raw difference (radians)
phase_shift = rad2deg(mod(phase_diff + pi, 2*pi) - pi); % in degrees, normalize to [-pi, pi]


%% --- Display result ---
fprintf("Fundamental frequency: %.2f Hz\n", f(k_fund));
fprintf("Phase(V) = %.4f rad (%.2f°)\n", phase_v, rad2deg(phase_v));
fprintf("Phase(I) = %.4f rad (%.2f°)\n", phase_i, rad2deg(phase_i));
fprintf("Phase Shift (V - I) = %.4f rad (%.4f°)\n", ...
        phase_v - phase_i, phase_shift);

%% --- Plot (optional) ---
figure;
subplot(2,1,1);
plot(t, v, 'b', t, i, 'r');
legend('Voltage','Current');
xlabel('Time (s)');
ylabel('Amplitude');
title('Voltage and Current Signals');
grid on;

subplot(2,1,2);
stem(f, V_mag, 'b'); hold on;
stem(f, I_mag, 'r');
xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Voltage','Current');
title('FFT Spectrum');
grid on;
