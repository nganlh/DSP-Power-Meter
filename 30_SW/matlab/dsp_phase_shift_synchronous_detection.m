clc; clear; close all;

fs = 7200;
f = 50;
N = round(fs / f * 3);
t = (0:N-1)/fs;
phase_deg = 30.1234;

% Example signals
v = 311*sin(2*pi*f*t);
i_pure = 311*sin(2*pi*f*t - deg2rad(phase_deg));

% Generate PWM waveform & apply to current
f_pwm = 20e3;                 % PWM frequency
pwm_carrier = square(2*pi*f_pwm*t, 50);   % 50% duty PWM: +1 / -1
pwm_carrier = (pwm_carrier + 1) / 2;  % Convert from [-1,1] → [0,1]
i = i_pure .* pwm_carrier;   % Modulated current

% References
ref_sin = sin(2*pi*f*t);
ref_cos = cos(2*pi*f*t);

% Voltage demodulation
Iv = (2/N)*sum(v .* ref_sin);
Qv = (2/N)*sum(v .* ref_cos);
phi_v = atan2(Qv, Iv);

% Current demodulation
Ii = (2/N)*sum(i .* ref_sin);
Qi = (2/N)*sum(i .* ref_cos);
phi_i = atan2(Qi, Ii);

% Phase shift
phase_shift = mod(phi_v - phi_i + pi, 2*pi) - pi;
fprintf('Phase shift = %.4f deg (expected: %.4f)\n', ...
    rad2deg(phase_shift), mod(phase_deg + 180, 360) - 180);



% === Plot ===
figure('Position',[100 100 1000 600]);

subplot(3,1,1);
plot(t*1000, v, 'b', 'LineWidth', 1.2);
title('Voltage Signal');
xlabel('Time (ms)'); ylabel('Amplitude (V)');
grid on;

subplot(3,1,2);
plot(t*1000, i_pure, 'r', 'LineWidth', 1.2);
title('Pure Current Signal');
xlabel('Time (ms)'); ylabel('Amplitude (A)');
grid on;

subplot(3,1,3);
plot(t*1000, i, 'k', 'LineWidth', 1);
hold on;
plot(t*1000, pwm_carrier * max(i_pure), 'g--');
title('PWM-Modulated Current');
xlabel('Time (ms)'); ylabel('Amplitude (A)');
legend('PWM Current', 'PWM Carrier (scaled)');
grid on;