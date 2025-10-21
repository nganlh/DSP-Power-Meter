clc; clear; close all;

% --- Parameters ---
fs = 7200;              % Sampling frequency (Hz)
f_signal = 50;          % Fundamental frequency (Hz)
N = 720;                % Number of samples (not power of 2 needed)
t = (0:N-1)/fs;         % Time vector

V_amp = 220*sqrt(2);    % Voltage amplitude (Vrms = 220V)
I_amp = 10*sqrt(2);     % Current amplitude (Irms = 10A)
phi_deg = -0.1;           % Current phase lag (deg)
phi = deg2rad(phi_deg); % Convert to radians

% --- Generate test signals ---
v = V_amp * sin(2*pi*f_signal*t);              % Voltage
i = I_amp * sin(2*pi*f_signal*t - phi);        % Current (lagging)

% Optional: Add offset or noise
v = v + 0.5*randn(size(v));
i = i + 0.5*randn(size(i));

% --- Remove DC offset ---
v = v - mean(v);
i = i - mean(i);

% --- Goertzel algorithm setup ---
k = round(f_signal * N / fs);      % Target frequency bin
omega = 2*pi*k / N;
coeff = 2*cos(omega);

% --- Goertzel for voltage ---
s_prev = 0; s_prev2 = 0;
for n = 1:N
    s = v(n) + coeff*s_prev - s_prev2;
    s_prev2 = s_prev;
    s_prev = s;
end
v_real = s_prev - s_prev2*cos(omega);
v_imag = s_prev2*sin(omega);
V_phasor = v_real + 1j*v_imag;

% --- Goertzel for current ---
s_prev = 0; s_prev2 = 0;
for n = 1:N
    s = i(n) + coeff*s_prev - s_prev2;
    s_prev2 = s_prev;
    s_prev = s;
end
i_real = s_prev - s_prev2*cos(omega);
i_imag = s_prev2*sin(omega);
I_phasor = i_real + 1j*i_imag;

% --- Phase calculation ---
phase_v = angle(V_phasor);
phase_i = angle(I_phasor);
phase_shift = wrapToPi(phase_v - phase_i);   % radians
phase_shift_deg = rad2deg(phase_shift);      % degrees

fprintf('Estimated phase shift = %.3f° (expected = %.1f°)\n', ...
    phase_shift_deg, phi_deg);

% --- Plot signals ---
subplot(2,1,1);
plot(t*1000, v, 'b', t*1000, i, 'r');
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Voltage','Current');
title('Input Signals');

subplot(2,1,2);
stem([real(V_phasor), imag(V_phasor)], 'b', 'filled'); hold on;
stem([real(I_phasor), imag(I_phasor)], 'r', 'filled');
legend('V phasor','I phasor');
xlabel('Real / Imag index'); ylabel('Magnitude');
title(sprintf('Goertzel Phase Shift: %.3f°', phase_shift_deg));
grid on;
