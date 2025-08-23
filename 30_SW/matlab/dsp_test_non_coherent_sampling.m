%% RMS & Leakage Demo
clear; clc; close all;

% --- Parameters ---
fs = 5000;            % sampling frequency (Hz)
f0 = 50;              % signal frequency (Hz)
Vrms_true = 220;      % target RMS value for comparison
A = Vrms_true * sqrt(2);  % peak amplitude for sine (A = Vrms*sqrt(2))

%% CASE A: Non-coherent sampling (non-integer cycles)
T_noncoh = 0.123;            % duration not integer number of cycles (intentionally)
tA = (0:1/fs:T_noncoh-1/fs)'; % column vector
xA = A * sin(2*pi*f0*tA);

% direct time-domain RMS (should be very close to true if sampling long enough)
rms_time_A = sqrt(mean(xA.^2));

% FFT amplitude method (single-sided amplitude, naive scaling)
N = length(xA);
X = fft(xA)/N;
amps = 2*abs(X(1:floor(N/2)));  % single-sided amplitude spectrum (except DC)
freqs = (0:floor(N/2)-1)*(fs/N);
% find fundamental bin
[~,kfund] = min(abs(freqs - f0));
amp_fund_naive = amps(kfund);   % naive amplitude from FFT
rms_fund_naive = amp_fund_naive / sqrt(2);

%% CASE B: Non-coherent + Hann window (no gain correction)
w = hann(N);          % column vector if N-by-1 signal; match dims
xB = xA .* w;
% FFT after window
Xw = fft(xB)/N;
amps_w = 2*abs(Xw(1:floor(N/2)));
amp_fund_w_naive = amps_w(kfund);
rms_fund_w_naive = amp_fund_w_naive / sqrt(2);

% Correct for coherent gain of hann window
win_coherent_gain = sum(w)/N;   % coherent gain = mean(window)
amp_fund_w_corrected = amp_fund_w_naive / win_coherent_gain;
rms_fund_w_corrected = amp_fund_w_corrected / sqrt(2);

%% CASE C: Coherent sampling (integer cycles) - choose T so M cycles
M = 10; % integer number of cycles
T_coh = M / f0; 
tC = (0:1/fs:T_coh-1/fs)'; 
xC = A * sin(2*pi*f0*tC);
N_c = length(xC);
% FFT amplitude (naive) for coherent case
Xc = fft(xC)/N_c;
amps_c = 2*abs(Xc(1:floor(N_c/2)));
freqs_c = (0:floor(N_c/2)-1)*(fs/N_c);
[~,kfund_c] = min(abs(freqs_c - f0));
amp_fund_c = amps_c(kfund_c);
rms_fund_c = amp_fund_c / sqrt(2);

% direct time-domain RMS for coherent
rms_time_C = sqrt(mean(xC.^2));

%% CASE D: Synchronous detection (correlation) on non-coherent case
% Using xA sample set (non-coherent). compute amplitude by correlation with
% sin/cos at exact f0 frequency (no leakage).
sin_ref = sin(2*pi*f0*tA);
cos_ref = cos(2*pi*f0*tA);
% In-phase and quadrature coefficients (peak amplitude scaling)
Xc_inphase = (2/N) * sum(xA .* cos_ref);
Xc_quad     = (2/N) * sum(xA .* sin_ref);
amp_sync = sqrt(Xc_inphase^2 + Xc_quad^2);   % peak amplitude
rms_sync = amp_sync / sqrt(2);

%% PRINT RESULTS
fprintf('True Vrms (design) = %.6f V\n\n', Vrms_true);

fprintf('CASE A (non-coherent, time-domain)    : Vrms_time = %.6f V   (should be close)\n', rms_time_A);
fprintf('CASE A (non-coherent, FFT naive)      : Vrms_fund = %.6f V   (underestimated due to leakage)\n\n', rms_fund_naive);

fprintf('CASE B (non-coherent + Hann, no corr) : Vrms_fund = %.6f V   (window attenuates amplitude)\n', rms_fund_w_naive);
fprintf('CASE B (Hann corrected by gain)       : Vrms_fund = %.6f V   (after dividing by window gain)\n\n', rms_fund_w_corrected);

fprintf('CASE C (coherent, FFT naive)          : Vrms_fund = %.6f V   (should match true)\n', rms_fund_c);
fprintf('CASE C (coherent, time-domain rms)    : Vrms_time  = %.6f V\n\n', rms_time_C);

fprintf('CASE D (synchronous detection)        : Vrms_sync  = %.6f V   (robust, little leakage effect)\n\n', rms_sync);

%% PLOTS
figure('Name','Signals and Spectra','Units','normalized','Position',[0.05 0.05 0.9 0.8]);

subplot(3,2,1);
plot(tA, xA); title('Case A: Non-coherent (time)'); xlabel('s'); ylabel('V'); xlim([0 min(0.02,tA(end))]);

subplot(3,2,2);
plot(freqs, amps); title('Case A: Spectrum (naive)'); xlabel('Hz'); ylabel('Amplitude'); xlim([0 500]);

subplot(3,2,3);
plot(tA, xB); title('Case B: Non-coherent * Hann'); xlabel('s'); ylabel('V'); xlim([0 min(0.02,tA(end))]);

subplot(3,2,4);
plot(freqs, amps_w); title('Case B: Spectrum (windowed)'); xlabel('Hz'); ylabel('Amplitude'); xlim([0 500]);

subplot(3,2,5);
plot(tC, xC); title('Case C: Coherent (integer cycles)'); xlabel('s'); ylabel('V'); xlim([0 min(0.02,tC(end))]);

subplot(3,2,6);
plot(freqs_c, amps_c); title('Case C: Spectrum (coherent)'); xlabel('Hz'); ylabel('Amplitude'); xlim([0 500]);

% Add an overall annotation instead of sgtitle for compatibility
annotation('textbox',[0.02 0.95 0.96 0.04],'String','RMS & Leakage Demonstration','EdgeColor','none','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
