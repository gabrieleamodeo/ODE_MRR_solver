%% Optical Ordinary Differential Equation (ODE) Solver via Silicon Photonic MRR
% Based on: "Compact tunable silicon photonic ODE solver", Optics Express (2014)
% Implementation of the Through-Port First-Order All-Optical ODE Solver
%
% This script compares:
%  1. Ideal First-Order ODE Transfer Function (Coupled-Mode Theory / Baseband Target ODE)
%  2. Full Periodic Physical MRR Transfer Function (Exact Add-Drop Microring Model)

clear; close all; clc;

%% 1. Physical Parameters
c          = 2.99792458e8;       % Speed of light in vacuum [m/s]
lambda0    = 1550.391e-9;        % Resonant wavelength [m] (1550.391 nm)
f0         = c / lambda0;        % Optical carrier frequency [Hz]
w0         = 2 * pi * f0;        % Optical carrier angular frequency [rad/s]
L          = 178.98e-6;          % Ring circumference [m] (178.98 um)
ng         = 4.1850;             % Waveguide group index
alpha_dB_cm= 8.0;                % Waveguide power loss factor [dB/cm]

% Convert waveguide power loss to round-trip power loss factor eta
alpha_dB_m = alpha_dB_cm * 100;                 % [dB/m]
alpha_lin  = alpha_dB_m * log(10) / 10;         % Linear power attenuation coefficient [1/m]
eta        = 1 - exp(-alpha_lin * L);           % Round-trip intrinsic power loss
a_loss     = sqrt(1 - eta);                     % Round-trip field amplitude transmission

% Initial power coupling coefficients of the two couplers
kappa1     = 0.08;                              % Coupler 1 (bus/through port)
kappa2     = 0.04;                              % Coupler 2 (drop port)
r1         = sqrt(1 - kappa1);                  % Coupler 1 field transmission coefficient
r2         = sqrt(1 - kappa2);                  % Coupler 2 field transmission coefficient

% Ring cavity round-trip time and Free Spectral Range (FSR)
tau_rt     = ng * L / c;                        % Round-trip time [s] (~2.498 ps)
FSR        = 1 / tau_rt;                        % Free Spectral Range [Hz] (~400.24 GHz)

%% 2. Coupled-Mode Theory (CMT) Formulation & ODE Coefficients
% Quality factors from CMT:
% Qi  = - (w0 * ng * L) / (c * ln(1 - eta))
% Qe1 = - (w0 * ng * L) / (c * ln(1 - kappa1))
% Qe2 = - (w0 * ng * L) / (c * ln(1 - kappa2))
Qi         = - (w0 * ng * L) / (c * log(1 - eta));
Qe1        = - (w0 * ng * L) / (c * log(1 - kappa1));
Qe2        = - (w0 * ng * L) / (c * log(1 - kappa2));

% Decay rates (gamma = w0 / (2*Q)):
gamma_i    = w0 / (2 * Qi);      % Intrinsic cavity decay rate [rad/s]
gamma_e1   = w0 / (2 * Qe1);     % Coupling decay rate to waveguide 1 [rad/s]
gamma_e2   = w0 / (2 * Qe2);     % Coupling decay rate to waveguide 2 [rad/s]

% Total quality factor and loaded 3-dB angular bandwidth
Q_total    = 1 / (1/Qi + 1/Qe1 + 1/Qe2);
Delta_w3dB = w0 / Q_total;       % [rad/s]

% ODE Constant Coefficients:
% Target ODE: dy(t)/dt + a0*y(t) = dx(t)/dt + b0*x(t)
a0         = gamma_i + gamma_e1 + gamma_e2;   % a0 = w0 / (2*Q_total) = Delta_w3dB / 2
b0         = gamma_i + gamma_e2 - gamma_e1;

fprintf('========================================================================\n');
fprintf('                Silicon Photonic ODE Solver Parameters                  \n');
fprintf('========================================================================\n');
fprintf('Resonant Wavelength (lambda0)       : %.3f nm\n', lambda0 * 1e9);
fprintf('Optical Carrier Frequency (f0)      : %.3f THz\n', f0 / 1e12);
fprintf('Ring Circumference (L)              : %.2f um\n', L * 1e6);
fprintf('Waveguide Group Index (ng)          : %.4f\n', ng);
fprintf('Round-trip Time (tau_rt)            : %.3f ps\n', tau_rt * 1e12);
fprintf('Free Spectral Range (FSR)           : %.2f GHz (%.3f nm)\n', FSR / 1e9, (lambda0^2 / (ng*L)) * 1e9);
fprintf('Waveguide Loss                      : %.2f dB/cm (Round-trip loss eta = %.4f)\n', alpha_dB_cm, eta);
fprintf('Coupling Coefficients               : kappa1 = %.4f, kappa2 = %.4f\n', kappa1, kappa2);
fprintf('------------------------------------------------------------------------\n');
fprintf('Cavity Quality Factors              : Qi = %.2e, Qe1 = %.2e, Qe2 = %.2e\n', Qi, Qe1, Qe2);
fprintf('Loaded Quality Factor (Q)           : %.2e\n', Q_total);
fprintf('Decay Rates (gamma)                 : gamma_i  = %.3e rad/s (%.2f GHz)\n', gamma_i, gamma_i/(2*pi*1e9));
fprintf('                                      gamma_e1 = %.3e rad/s (%.2f GHz)\n', gamma_e1, gamma_e1/(2*pi*1e9));
fprintf('                                      gamma_e2 = %.3e rad/s (%.2f GHz)\n', gamma_e2, gamma_e2/(2*pi*1e9));
fprintf('ODE Coefficient a0                  : %.4e rad/s (Bandwidth: %.2f GHz)\n', a0, a0/(pi*1e9));
fprintf('ODE Coefficient b0                  : %.4e rad/s\n', b0);
if b0 < 0
    fprintf('Coupling Regime                     : Over-coupled (b0 < 0)\n');
elseif b0 == 0
    fprintf('Coupling Regime                     : Critically coupled (b0 = 0, Differentiator)\n');
else
    fprintf('Coupling Regime                     : Under-coupled (b0 > 0)\n');
end
fprintf('Through-port Static Notch (b0/a0)^2 : %.4f (%.2f dB)\n', (b0/a0)^2, 10*log10((b0/a0)^2));
fprintf('========================================================================\n\n');

%% 3. Signal Definition: 10 Gb/s Gaussian-like Optical Pulse (FWHM = 45 ps)
FWHM       = 45e-12;             % Full Width at Half Maximum: 45 ps (~40-50 ps)
t0         = 0.3e-9;             % Center of the pulse [s] (300 ps)

% Simulation time grid
t_start    = -0.5e-9;            % -0.5 ns
t_end      = 1.5e-9;             % 1.5 ns
N          = 65536;              % Number of grid points (power of 2 for FFT)
t          = linspace(t_start, t_end, N);
dt         = t(2) - t(1);

% Optical complex envelope x(t)
% Intensity: I_in(t) = exp(-4*log(2)*((t - t0)/FWHM)^2) -> Field: x(t) = sqrt(I_in(t))
x_fun      = @(t_val) exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);
x          = x_fun(t);           % Input field envelope

%% 4. Frequency-Domain Optical Filtering (FFT/IFFT)
fs         = 1 / dt;                          % Sampling frequency [Hz]
f          = linspace(-fs/2, fs/2 - fs/N, N); % Centered baseband frequency axis [Hz]
w          = 2 * pi * f;                      % Baseband angular frequency [rad/s]

% (A) Ideal First-Order ODE Transfer Function (Target ODE via CMT)
%     T_ODE(w) = (j*w + b0) / (j*w + a0)
T_ode      = (1j * w + b0) ./ (1j * w + a0);

% (B) Full Periodic Physical MRR Through-Port Transfer Function
%     H_MRR(w) = (r1 - r2*a*exp(-j*w*tau_rt)) / (1 - r1*r2*a*exp(-j*w*tau_rt))
H_mrr      = (r1 - r2 * a_loss * exp(-1j * w * tau_rt)) ./ ...
             (1 - r1 * r2 * a_loss * exp(-1j * w * tau_rt));

% Optical filtering via FFT
X_w        = fftshift(fft(x));           % Input pulse spectrum
Y_ode_w    = X_w .* T_ode;               % Output spectrum via Ideal ODE filter
Y_mrr_w    = X_w .* H_mrr;               % Output spectrum via Full Physical MRR

y_ode      = ifft(ifftshift(Y_ode_w));   % Time-domain field via Ideal ODE filter
y_mrr      = ifft(ifftshift(Y_mrr_w));   % Time-domain field via Physical MRR

%% 5. Power Intensities & Normalized Mean Square Error (NMSE)
I_in       = abs(x).^2;                  % Input intensity
I_ode      = abs(y_ode).^2;              % Ideal ODE output intensity
I_mrr      = abs(y_mrr).^2;              % Physical MRR solver output intensity

% Normalized Mean Square Error (NMSE) between Ideal ODE and Physical MRR
NMSE       = sum((I_ode - I_mrr).^2) / sum(I_ode.^2);
NMSE_dB    = 10 * log10(NMSE);

fprintf('========================================================================\n');
fprintf('                         Verification Results                           \n');
fprintf('========================================================================\n');
fprintf('Signal Type                         : 10 Gb/s Gaussian Pulse (FWHM = 45 ps)\n');
fprintf('Pulse FWHM                          : %.2f ps\n', FWHM * 1e12);
fprintf('Normalized Mean Square Error (NMSE) : %.4e\n', NMSE);
fprintf('NMSE in Percentage                  : %.6f %%\n', NMSE * 100);
fprintf('NMSE in dB                          : %.2f dB\n', NMSE_dB);
fprintf('========================================================================\n\n');

%% 6. Visualization & Plotting
set(0, 'DefaultAxesFontSize', 11, 'DefaultLineLineWidth', 1.8);

% -------------------------------------------------------------------------
% Figure 1: Time-Domain Waveforms and Physical MRR Verification
% -------------------------------------------------------------------------
figure('Name', 'Photonic ODE Solver - Time Domain Comparison', 'Position', [80, 80, 950, 820], 'Color', 'w');

% Subplot 1: Input Waveform
subplot(3, 1, 1);
plot(t * 1e12, I_in, 'Color', [0.1, 0.1, 0.1], 'DisplayName', 'Input Intensity |x(t)|^2');
hold on;
plot(t * 1e12, x, '--', 'Color', [0.3, 0.5, 0.9], 'LineWidth', 1.4, 'DisplayName', 'Input Field Envelope x(t)');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Amplitude / Intensity');
title('(a) Input 10 Gb/s Gaussian Pulse (FWHM = 45 ps)');
legend('Location', 'northeast');

% Subplot 2: Output Intensity Comparison
subplot(3, 1, 2);
plot(t * 1e12, I_ode, 'b-', 'DisplayName', 'Ideal ODE Target: |y_{ODE}(t)|^2');
hold on;
plot(t * 1e12, I_mrr, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Full Physical MRR: |y_{MRR}(t)|^2');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Output Intensity [a.u.]');
title(sprintf('(b) Output Intensity Comparison: Ideal Target ODE vs. Physical MRR (NMSE = %.2e)', NMSE));
legend('Location', 'northeast');

% Subplot 3: Physical Error & Residuals
subplot(3, 1, 3);
plot(t * 1e12, abs(I_ode - I_mrr), 'r-', 'DisplayName', 'Physical Deviation: |I_{ODE}(t) - I_{MRR}(t)|');
hold on;
plot(t * 1e12, real(y_ode), 'b:', 'LineWidth', 1.4, 'DisplayName', 'Re\{y_{ODE}(t)\}');
plot(t * 1e12, real(y_mrr), 'm-.', 'LineWidth', 1.2, 'DisplayName', 'Re\{y_{MRR}(t)\}');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Field / Intensity Error');
title('(c) Physical Deviation & Real Field Envelopes');
legend('Location', 'northeast');

% -------------------------------------------------------------------------
% Figure 2: Baseband Frequency Response Comparison (Zoom: +/- 50 GHz)
% -------------------------------------------------------------------------
figure('Name', 'Photonic ODE Solver - Baseband Frequency Response', 'Position', [120, 120, 950, 620], 'Color', 'w');

f_GHz = f / 1e9;
freq_limit = 50; % Plot range +/- 50 GHz

subplot(2, 1, 1);
plot(f_GHz, 10 * log10(abs(T_ode).^2), 'b--', 'LineWidth', 2.0, 'DisplayName', 'Ideal Target ODE |T_{ODE}(\omega)|^2');
hold on;
plot(f_GHz, 10 * log10(abs(H_mrr).^2), 'r-', 'LineWidth', 1.6, 'DisplayName', 'Full Physical MRR |H_{MRR}(\omega)|^2');
input_spec_dB = 10 * log10(abs(X_w).^2 / max(abs(X_w).^2));
plot(f_GHz, input_spec_dB, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Input Pulse Spectrum (Normalized)');
yline(10 * log10((b0/a0)^2), 'g-.', 'LineWidth', 1.2, 'DisplayName', sprintf('Resonance Notch (|b_0/a_0|^2 = %.2f dB)', 10*log10((b0/a0)^2)));
grid on; box on;
xlim([-freq_limit, freq_limit]);
ylim([-35, 5]);
xlabel('Frequency Offset f - f_0 [GHz]');
ylabel('Magnitude [dB]');
title('(a) Baseband Magnitude Response: Ideal ODE vs. Full Physical MRR');
legend('Location', 'southeast');

subplot(2, 1, 2);
plot(f_GHz, angle(T_ode) * (180 / pi), 'b--', 'LineWidth', 2.0, 'DisplayName', '\angle T_{ODE}(\omega) (Ideal ODE)');
hold on;
plot(f_GHz, angle(H_mrr) * (180 / pi), 'r-', 'LineWidth', 1.6, 'DisplayName', '\angle H_{MRR}(\omega) (Physical MRR)');
grid on; box on;
xlim([-freq_limit, freq_limit]);
ylim([-190, 190]);
xlabel('Frequency Offset f - f_0 [GHz]');
ylabel('Phase [degrees]');
title('(b) Baseband Phase Response Comparison');
legend('Location', 'southeast');

% -------------------------------------------------------------------------
% Figure 3: Broadband Spectral View & Free Spectral Range (FSR) Periodicity
% -------------------------------------------------------------------------
figure('Name', 'Photonic ODE Solver - Broadband Spectrum & FSR Periodicity', 'Position', [160, 160, 950, 520], 'Color', 'w');

broad_limit = 600; % +/- 600 GHz (showing > 1 FSR)

plot(f_GHz, 10 * log10(abs(H_mrr).^2), 'r-', 'LineWidth', 1.8, 'DisplayName', 'Full Periodic Physical MRR |H_{MRR}(\omega)|^2');
hold on;
plot(f_GHz, 10 * log10(abs(T_ode).^2), 'b--', 'LineWidth', 1.6, 'DisplayName', 'Single-Resonance Ideal ODE |T_{ODE}(\omega)|^2');
plot(f_GHz, input_spec_dB, 'k:', 'LineWidth', 1.3, 'DisplayName', 'Input Pulse Spectrum (10 Gb/s)');
xline(-FSR/1e9, 'm-.', 'LineWidth', 1.2, 'DisplayName', sprintf('-FSR (-%.1f GHz)', FSR/1e9));
xline(+FSR/1e9, 'm-.', 'LineWidth', 1.2, 'DisplayName', sprintf('+FSR (+%.1f GHz)', FSR/1e9));
grid on; box on;
xlim([-broad_limit, broad_limit]);
ylim([-35, 5]);
xlabel('Frequency Offset f - f_0 [GHz]');
ylabel('Magnitude [dB]');
title(sprintf('Broadband Spectral View: Physical MRR FSR Periodicity (FSR = %.2f GHz) vs. Baseband ODE Model', FSR/1e9));
legend('Location', 'southeast');
