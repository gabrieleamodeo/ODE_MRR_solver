
clear; close all; clc;

%% 1. Baseline Physical Parameters (from Wu et al., Opt. Express 2014)
c           = 2.99792458e8;       % Speed of light in vacuum [m/s]
lambda0     = 1550.391e-9;        % Target resonant wavelength [m] (1550.391 nm)
f0          = c / lambda0;        % Optical carrier frequency [Hz]
w0          = 2 * pi * f0;        % Optical carrier angular frequency [rad/s]
L           = 178.98e-6;          % Ring circumference [m] (178.98 um)
ng          = 4.1850;             % Waveguide group index
alpha_dB_cm = 8.0;                % Waveguide power loss factor [dB/cm]

% Waveguide loss conversions
alpha_dB_m  = alpha_dB_cm * 100;                 % [dB/m]
alpha_lin   = alpha_dB_m * log(10) / 10;         % Linear power attenuation coefficient [1/m]
eta         = 1 - exp(-alpha_lin * L);           % Round-trip intrinsic power loss
a_loss      = sqrt(1 - eta);                     % Round-trip field amplitude transmission

% Ring round-trip time and FSR
tau_rt      = ng * L / c;                        % Round-trip time [s] (~2.498 ps)
FSR         = 1 / tau_rt;                        % Free Spectral Range [Hz] (~400.24 GHz)

%% 2. Thermo-Optic Model & Coupler Phase Shifters
% Silicon thermo-optic coefficient around 1550 nm:
% dn/dT approx 1.86e-4 K^-1 (TE fundamental mode in SOI waveguide)
dn_dT       = 1.86e-4;                           % [1/K]

% Directional coupler splitting ratio at MZI coupler branches (from paper Section 3)
kappa0      = 0.0441;                            % Directional coupler power coupling coefficient
kappa_max   = 4 * kappa0 * (1 - kappa0);         % Max achievable MZI power coupling (~0.1686)

% Physical lengths of bus arm microheaters (Wu et al. 2014):
% Heater 1: Lb1 = 116.80 um; Heater 2: Lb2 = 47.12 um
Lb1         = 116.80e-6;                         % Bus arm 1 heater length [m]
Lb2         = 47.12e-6;                          % Bus arm 2 heater length [m]

% Target coupling coefficients for the target ODE:
% dy(t)/dt + a0*y(t) = dx(t)/dt + b0*x(t)
kappa1_target = 0.08;                            % Bus/through port coupler
kappa2_target = 0.04;                            % Drop port coupler

% Analytical inversion of MZI coupling formula:
% kappa = 2*kappa0*(1-kappa0)*(1 + cos(Delta_phi))
cos_dphi1   = (kappa1_target / (2 * kappa0 * (1 - kappa0))) - 1;
cos_dphi2   = (kappa2_target / (2 * kappa0 * (1 - kappa0))) - 1;
dphi1_nom   = acos(min(max(cos_dphi1, -1), 1));  % Nominal phase difference coupler 1 [rad]
dphi2_nom   = acos(min(max(cos_dphi2, -1), 1));  % Nominal phase difference coupler 2 [rad]

%% 3. Thermal Noise / Error Configuration
% Standard laboratory temperature control precision (Standard Deviation sigma_T):
% You can change sigma_T to match your lab equipment!
% Typical values: 0.005 (5 mK, high-end PID), 0.020 (20 mK, standard), 0.050 (50 mK, open loop)
sigma_T_heater = 0.020;  % [K] Standard deviation of heater temperature control (20 mK)
sigma_T_ring   = 0.015;  % [K] Standard deviation of ring substrate / cross-talk drift (15 mK)
thermal_crosstalk = 0.06; % 6% thermal leakage from heaters to the microring cavity

fprintf('========================================================================\n');
fprintf('       Photonic MRR ODE Solver - Thermal Noise & Control Error          \n');
fprintf('========================================================================\n');
fprintf('Carrier Wavelength (lambda0)        : %.3f nm\n', lambda0 * 1e9);
fprintf('Waveguide Thermo-Optic (dn/dT)      : %.2e K^-1\n', dn_dT);
fprintf('Nominal Couplers                    : kappa1 = %.4f, kappa2 = %.4f\n', kappa1_target, kappa2_target);
fprintf('Nominal Coupler Phase Shifts        : dphi1 = %.2f deg, dphi2 = %.2f deg\n', dphi1_nom*180/pi, dphi2_nom*180/pi);
fprintf('Heater 1 Length (Lb1)               : %.2f um\n', Lb1 * 1e6);
fprintf('Heater 2 Length (Lb2)               : %.2f um\n', Lb2 * 1e6);
fprintf('Phase Shifter Thermal Noise (sigma) : %.2f mK (%.4f K)\n', sigma_T_heater * 1e3, sigma_T_heater);
fprintf('Cavity Resonance Drift (sigma)      : %.2f mK (%.4f K)\n', sigma_T_ring * 1e3, sigma_T_ring);
fprintf('========================================================================\n\n');

%% 4. Nominal ODE & MRR Transfer Functions (Ideal / Noiseless)
r1_nom      = sqrt(1 - kappa1_target);
r2_nom      = sqrt(1 - kappa2_target);

% Quality factors and decay rates (Nominal)
Qi_nom      = - (w0 * ng * L) / (c * log(1 - eta));
Qe1_nom     = - (w0 * ng * L) / (c * log(1 - kappa1_target));
Qe2_nom     = - (w0 * ng * L) / (c * log(1 - kappa2_target));
gamma_i     = w0 / (2 * Qi_nom);
gamma_e1_nom= w0 / (2 * Qe1_nom);
gamma_e2_nom= w0 / (2 * Qe2_nom);

a0_target   = gamma_i + gamma_e1_nom + gamma_e2_nom;
b0_target   = gamma_i + gamma_e2_nom - gamma_e1_nom;

%% 5. Signal Definition: 10 Gb/s Gaussian Pulse (FWHM = 45 ps)
FWHM        = 45e-12;             % 45 ps
t0          = 0.3e-9;             % Center: 300 ps
t_start     = -0.5e-9;
t_end       = 1.5e-9;
N           = 65536;              % Power of 2 grid
t           = linspace(t_start, t_end, N);
dt          = t(2) - t(1);
fs          = 1 / dt;
f           = linspace(-fs/2, fs/2 - fs/N, N); % Centered baseband frequency [Hz]
w           = 2 * pi * f;

% Gaussian input pulse
x_fun       = @(t_val) exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);
x           = x_fun(t);
I_in        = abs(x).^2;
X_w         = fftshift(fft(x));

% Ideal Target ODE Transfer Function: (j*w + b0) / (j*w + a0)
T_ode_ideal = (1j * w + b0_target) ./ (1j * w + a0_target);
Y_ode_w     = X_w .* T_ode_ideal;
y_ode_ideal = ifft(ifftshift(Y_ode_w));
I_ode_ideal = abs(y_ode_ideal).^2;

% Nominal MRR (Noiseless)
H_mrr_nom   = (r1_nom - r2_nom * a_loss * exp(-1j * w * tau_rt)) ./ ...
              (1 - r1_nom * r2_nom * a_loss * exp(-1j * w * tau_rt));
Y_mrr_nom_w = X_w .* H_mrr_nom;
y_mrr_nom   = ifft(ifftshift(Y_mrr_nom_w));
I_mrr_nom   = abs(y_mrr_nom).^2;

NMSE_nominal = sum((I_ode_ideal - I_mrr_nom).^2) / sum(I_ode_ideal.^2);

%% 6. Single Perturbed Realization (Realistic Lab Shot)
rng(42); % Fixed seed for reproducible single realization display

% Random temperature errors
dT_h1       = randn * sigma_T_heater;
dT_h2       = randn * sigma_T_heater;
dT_cavity   = randn * sigma_T_ring + thermal_crosstalk * (dT_h1 + dT_h2) / 2;

% Thermo-optic phase perturbations:
% delta_phi = (2*pi / lambda0) * dn/dT * L_heater * dT
dphi1_noisy = dphi1_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dT_h1;
dphi2_noisy = dphi2_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dT_h2;

% Perturbed power coupling coefficients
kappa1_noisy= 2 * kappa0 * (1 - kappa0) * (1 + cos(dphi1_noisy));
kappa2_noisy= 2 * kappa0 * (1 - kappa0) * (1 + cos(dphi2_noisy));
% Enforce physical bounds [0, 1)
kappa1_noisy= min(max(kappa1_noisy, 1e-6), 0.99);
kappa2_noisy= min(max(kappa2_noisy, 1e-6), 0.99);
r1_noisy    = sqrt(1 - kappa1_noisy);
r2_noisy    = sqrt(1 - kappa2_noisy);

% Microring round-trip phase detuning:
% delta_phi_rt = (2*pi / lambda0) * dn/dT * L * dT_cavity
dphi_rt_noisy = (2 * pi / lambda0) * dn_dT * L * dT_cavity;
% Corresponding resonance frequency shift: Delta_f = dphi_rt / (2*pi*tau_rt)
Delta_f_res   = dphi_rt_noisy / (2 * pi * tau_rt);

% Perturbed Physical MRR Transfer Function
% Round-trip phase factor: exp(-j * (w * tau_rt + dphi_rt_noisy))
H_mrr_noisy = (r1_noisy - r2_noisy * a_loss * exp(-1j * (w * tau_rt + dphi_rt_noisy))) ./ ...
              (1 - r1_noisy * r2_noisy * a_loss * exp(-1j * (w * tau_rt + dphi_rt_noisy)));

% Output spectrum and temporal intensity
Y_mrr_noisy_w = X_w .* H_mrr_noisy;
y_mrr_noisy   = ifft(ifftshift(Y_mrr_noisy_w));
I_mrr_noisy   = abs(y_mrr_noisy).^2;

% Realized ODE coefficients from noisy Q factors
Qe1_noisy   = - (w0 * ng * L) / (c * log(1 - kappa1_noisy));
Qe2_noisy   = - (w0 * ng * L) / (c * log(1 - kappa2_noisy));
gamma_e1_n  = w0 / (2 * Qe1_noisy);
gamma_e2_n  = w0 / (2 * Qe2_noisy);
a0_noisy    = gamma_i + gamma_e1_n + gamma_e2_n;
b0_noisy    = gamma_i + gamma_e2_n - gamma_e1_n;

NMSE_noisy  = sum((I_ode_ideal - I_mrr_noisy).^2) / sum(I_ode_ideal.^2);

fprintf('--- Single Realization Error Breakdown ---\n');
fprintf('Heater 1 Error (dT1)                : %+.2f mK -> kappa1 = %.4f (nom: %.4f, delta: %+.3e)\n', ...
    dT_h1 * 1e3, kappa1_noisy, kappa1_target, kappa1_noisy - kappa1_target);
fprintf('Heater 2 Error (dT2)                : %+.2f mK -> kappa2 = %.4f (nom: %.4f, delta: %+.3e)\n', ...
    dT_h2 * 1e3, kappa2_noisy, kappa2_target, kappa2_noisy - kappa2_target);
fprintf('Cavity Temp Drift (dT_cav)          : %+.2f mK -> Resonance Shift: %+.2f MHz\n', ...
    dT_cavity * 1e3, Delta_f_res / 1e6);
fprintf('Target ODE Coeffs                   : a0 = %.4e rad/s, b0 = %.4e rad/s\n', a0_target, b0_target);
fprintf('Realized Noisy ODE Coeffs           : a0 = %.4e rad/s, b0 = %.4e rad/s\n', a0_noisy, b0_noisy);
fprintf('Nominal MRR NMSE (Intrinsic Model)  : %.4e (%.2f dB)\n', NMSE_nominal, 10*log10(NMSE_nominal));
fprintf('Thermally Perturbed MRR NMSE        : %.4e (%.2f dB)\n', NMSE_noisy, 10*log10(NMSE_noisy));
fprintf('NMSE Degradation due to Temp Error  : %+.2f dB\n', 10*log10(NMSE_noisy) - 10*log10(NMSE_nominal));
fprintf('========================================================================\n\n');

%% 7. Monte Carlo Statistical Simulation (N_trials Realizations)
N_trials     = 200;
NMSE_mc      = zeros(N_trials, 1);
I_mrr_mc     = zeros(N, N_trials);
a0_mc        = zeros(N_trials, 1);
b0_mc        = zeros(N_trials, 1);
Delta_f_mc   = zeros(N_trials, 1);

fprintf('Running Monte Carlo simulation (%d trials)...\n', N_trials);
for trial = 1:N_trials
    % Sample independent thermal variations
    dt1 = randn * sigma_T_heater;
    dt2 = randn * sigma_T_heater;
    dt_cav = randn * sigma_T_ring + thermal_crosstalk * (dt1 + dt2) / 2;

    % Phase & coupling
    dp1 = dphi1_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt1;
    dp2 = dphi2_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt2;
    k1  = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp1)), 1e-6), 0.99);
    k2  = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp2)), 1e-6), 0.99);
    r1_t = sqrt(1 - k1);
    r2_t = sqrt(1 - k2);

    % Cavity detuning
    dp_rt = (2 * pi / lambda0) * dn_dT * L * dt_cav;
    Delta_f_mc(trial) = dp_rt / (2 * pi * tau_rt);

    % Transfer function & pulse output
    H_t = (r1_t - r2_t * a_loss * exp(-1j * (w * tau_rt + dp_rt))) ./ ...
          (1 - r1_t * r2_t * a_loss * exp(-1j * (w * tau_rt + dp_rt)));
    y_t = ifft(ifftshift(X_w .* H_t));
    I_t = abs(y_t).^2;

    I_mrr_mc(:, trial) = I_t;
    NMSE_mc(trial)     = sum((I_ode_ideal - I_t).^2) / sum(I_ode_ideal.^2);

    % ODE coefficients
    qe1 = - (w0 * ng * L) / (c * log(1 - k1));
    qe2 = - (w0 * ng * L) / (c * log(1 - k2));
    a0_mc(trial) = gamma_i + (w0 / (2 * qe1)) + (w0 / (2 * qe2));
    b0_mc(trial) = gamma_i + (w0 / (2 * qe2)) - (w0 / (2 * qe1));
end

mean_NMSE    = mean(NMSE_mc);
std_NMSE     = std(NMSE_mc);
median_NMSE  = median(NMSE_mc);
mean_I_mc    = mean(I_mrr_mc, 2);
std_I_mc     = std(I_mrr_mc, 0, 2);

fprintf('Monte Carlo Complete:\n');
fprintf('  Mean NMSE   : %.4e (%.2f dB)\n', mean_NMSE, 10*log10(mean_NMSE));
fprintf('  Median NMSE : %.4e (%.2f dB)\n', median_NMSE, 10*log10(median_NMSE));
fprintf('  Worst NMSE  : %.4e (%.2f dB)\n', max(NMSE_mc), 10*log10(max(NMSE_mc)));
fprintf('  Best NMSE   : %.4e (%.2f dB)\n', min(NMSE_mc), 10*log10(min(NMSE_mc)));
fprintf('========================================================================\n\n');

%% 8. Sensitivity Analysis: NMSE vs. Thermal Control Accuracy (sigma_T)
fprintf('Running Sensitivity Sweep (NMSE vs. sigma_T)...\n');
sigma_T_vec  = [1, 2, 5, 10, 15, 20, 30, 40, 50, 75, 100] * 1e-3; % 1 mK to 100 mK
sweep_trials = 60;
NMSE_sweep_mean_dB = zeros(length(sigma_T_vec), 1);
NMSE_sweep_std_dB  = zeros(length(sigma_T_vec), 1);

for s_idx = 1:length(sigma_T_vec)
    sig_T = sigma_T_vec(s_idx);
    nmse_temp = zeros(sweep_trials, 1);
    for t_idx = 1:sweep_trials
        dt1 = randn * sig_T;
        dt2 = randn * sig_T;
        dt_cav = randn * (sig_T * 0.75) + thermal_crosstalk * (dt1 + dt2) / 2;

        dp1 = dphi1_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt1;
        dp2 = dphi2_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt2;
        k1  = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp1)), 1e-6), 0.99);
        k2  = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp2)), 1e-6), 0.99);
        r1_t = sqrt(1 - k1);
        r2_t = sqrt(1 - k2);

        dp_rt = (2 * pi / lambda0) * dn_dT * L * dt_cav;
        H_t = (r1_t - r2_t * a_loss * exp(-1j * (w * tau_rt + dp_rt))) ./ ...
              (1 - r1_t * r2_t * a_loss * exp(-1j * (w * tau_rt + dp_rt)));
        y_t = ifft(ifftshift(X_w .* H_t));
        nmse_temp(t_idx) = sum((I_ode_ideal - abs(y_t).^2).^2) / sum(I_ode_ideal.^2);
    end
    nmse_temp_dB = 10 * log10(nmse_temp);
    NMSE_sweep_mean_dB(s_idx) = mean(nmse_temp_dB);
    NMSE_sweep_std_dB(s_idx)  = std(nmse_temp_dB);
end

%% 9. Visualization & Plotting
set(0, 'DefaultAxesFontSize', 10, 'DefaultLineLineWidth', 1.6);

% -------------------------------------------------------------------------
% Figure 1: Single Realization Comparison (Time Domain)
% -------------------------------------------------------------------------
fig1 = figure('Name', 'Thermal Error - Time Domain Comparison', 'Position', [60, 60, 960, 820], 'Color', 'w');

subplot(3, 1, 1);
plot(t * 1e12, I_in, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Input Pulse |x(t)|^2');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]'); ylabel('Intensity');
title('(a) Input 10 Gb/s Gaussian-like Optical Pulse (FWHM = 45 ps)');
legend('Location', 'northeast');

subplot(3, 1, 2);
plot(t * 1e12, I_ode_ideal, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Ideal ODE Target');
hold on;
plot(t * 1e12, I_mrr_nom, 'g--', 'LineWidth', 1.6, 'DisplayName', sprintf('Nominal MRR (NMSE = %.1e)', NMSE_nominal));
plot(t * 1e12, I_mrr_noisy, 'r-.', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('MRR with Thermal Noise (\\sigma_T = 20 mK, NMSE = %.2e)', NMSE_noisy));
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]'); ylabel('Output Intensity [a.u.]');
title('(b) Time-Domain Output: Ideal vs. Nominal vs. Thermally Perturbed MRR');
legend('Location', 'northeast');

subplot(3, 1, 3);
plot(t * 1e12, abs(I_ode_ideal - I_mrr_nom), 'g--', 'LineWidth', 1.4, 'DisplayName', 'Error: |I_{ODE} - I_{MRR,nominal}|');
hold on;
plot(t * 1e12, abs(I_ode_ideal - I_mrr_noisy), 'r-', 'LineWidth', 1.6, 'DisplayName', 'Error: |I_{ODE} - I_{MRR,noisy}|');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]'); ylabel('Intensity Deviation');
title('(c) Instantaneous Calculation Deviation Due to Thermal Imprecision');
legend('Location', 'northeast');

% -------------------------------------------------------------------------
% Figure 2: Baseband Spectral Response & Detuning Effect (Zoom +/- 20 GHz)
% -------------------------------------------------------------------------
fig2 = figure('Name', 'Thermal Error - Baseband Spectral Response', 'Position', [100, 100, 960, 620], 'Color', 'w');
f_GHz = f / 1e9;
freq_limit = 20; % +/- 20 GHz (zoomed for clear notch visibility)

subplot(2, 1, 1);
plot(f_GHz, 10 * log10(abs(T_ode_ideal).^2), 'b--', 'LineWidth', 2.0, 'DisplayName', 'Ideal Target ODE');
hold on;
plot(f_GHz, 10 * log10(abs(H_mrr_nom).^2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Nominal MRR (Noiseless)');
plot(f_GHz, 10 * log10(abs(H_mrr_noisy).^2), 'r-', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Perturbed MRR (\\Delta f_{res} = %+.1f MHz)', Delta_f_res/1e6));
grid on; box on;
xlim([-freq_limit, freq_limit]); ylim([-35, 5]);
xlabel('Frequency Offset f - f_0 [GHz]'); ylabel('Transmission [dB]');
title('(a) Baseband Magnitude Response: Resonance Notch Shift and Asymmetry');
legend('Location', 'southeast');

subplot(2, 1, 2);
plot(f_GHz, angle(T_ode_ideal) * 180 / pi, 'b--', 'LineWidth', 2.0, 'DisplayName', 'Ideal Target ODE');
hold on;
plot(f_GHz, angle(H_mrr_nom) * 180 / pi, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Nominal MRR');
plot(f_GHz, angle(H_mrr_noisy) * 180 / pi, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Perturbed MRR');
grid on; box on;
xlim([-freq_limit, freq_limit]); ylim([-190, 190]);
xlabel('Frequency Offset f - f_0 [GHz]'); ylabel('Phase [deg]');
title('(b) Baseband Phase Response under Thermal Perturbation');
legend('Location', 'southeast');

% -------------------------------------------------------------------------
% Figure 3: Monte Carlo Statistical Analysis
% -------------------------------------------------------------------------
fig3 = figure('Name', 'Thermal Error - Monte Carlo Statistics', 'Position', [140, 140, 1000, 650], 'Color', 'w');

% Subplot 1: Waveform Uncertainty Bands
subplot(2, 2, [1, 3]);
t_ps = t * 1e12;
y_upper2 = (mean_I_mc + 2*std_I_mc)';
y_lower2 = (max(mean_I_mc - 2*std_I_mc, 0))';
y_upper1 = (mean_I_mc + std_I_mc)';
y_lower1 = (max(mean_I_mc - std_I_mc, 0))';

fill([t_ps, fliplr(t_ps)], [y_upper2, fliplr(y_lower2)], ...
    [1.0, 0.85, 0.85], 'EdgeColor', 'none', 'DisplayName', '\pm 2\sigma Band (95.4%)');
hold on;
fill([t_ps, fliplr(t_ps)], [y_upper1, fliplr(y_lower1)], ...
    [1.0, 0.65, 0.65], 'EdgeColor', 'none', 'DisplayName', '\pm 1\sigma Band (68.3%)');
plot(t_ps, I_ode_ideal, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Ideal ODE Target');
plot(t_ps, mean_I_mc, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Monte Carlo Mean Waveform');
grid on; box on;
xlim([100, 500]);
xlabel('Time [ps]'); ylabel('Output Intensity [a.u.]');
title(sprintf('(a) Waveform Uncertainty Bands (%d Monte Carlo Trials, \\sigma_T = 20 mK)', N_trials));
legend('Location', 'northeast');

% Subplot 2: NMSE Histogram
subplot(2, 2, 2);
histogram(10 * log10(NMSE_mc), 25, 'FaceColor', [0.2, 0.5, 0.8], 'EdgeColor', 'k');
xline(10 * log10(NMSE_nominal), 'g--', 'LineWidth', 2.0, 'DisplayName', sprintf('Nominal: %.1f dB', 10*log10(NMSE_nominal)));
xline(10 * log10(mean_NMSE), 'r-', 'LineWidth', 2.0, 'DisplayName', sprintf('Mean: %.1f dB', 10*log10(mean_NMSE)));
grid on; box on;
xlabel('NMSE [dB]'); ylabel('Occurrences');
title('(b) NMSE Distribution under Thermal Noise');
legend('Location', 'northwest');

% Subplot 3: ODE Parameter Scatter (a0 vs b0)
subplot(2, 2, 4);
scatter(a0_mc / 1e10, b0_mc / 1e10, 25, 10*log10(NMSE_mc), 'filled');
hold on;
plot(a0_target / 1e10, b0_target / 1e10, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'Target (a_0, b_0)');
colorbar; ylabel(colorbar, 'NMSE [dB]');
grid on; box on;
xlabel('ODE Coefficient a_0 [\times 10^{10} rad/s]');
ylabel('ODE Coefficient b_0 [\times 10^{10} rad/s]');
title('(c) Realized ODE Parameter Drift Cluster');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 4: Sensitivity Curve (NMSE vs Thermal Control Accuracy sigma_T)
% -------------------------------------------------------------------------
fig4 = figure('Name', 'Thermal Error - Lab Accuracy Sensitivity Curve', 'Position', [180, 180, 850, 520], 'Color', 'w');

sig_T_mK = sigma_T_vec * 1e3;
nmse_upper_dB = (NMSE_sweep_mean_dB + NMSE_sweep_std_dB)';
nmse_lower_dB = (NMSE_sweep_mean_dB - NMSE_sweep_std_dB)';

fill([sig_T_mK, fliplr(sig_T_mK)], [nmse_upper_dB, fliplr(nmse_lower_dB)], ...
    [0.85, 0.9, 1.0], 'EdgeColor', 'none', 'DisplayName', '\pm 1\sigma Spread');
hold on;
plot(sig_T_mK, NMSE_sweep_mean_dB, 'b-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'b', 'DisplayName', 'Mean NMSE [dB]');
yline(10 * log10(NMSE_nominal), 'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('Intrinsic Ring Approximation Limit (%.1f dB)', 10*log10(NMSE_nominal)));

% Mark typical laboratory benchmarks
xline(5,  'm-.', 'LineWidth', 1.4, 'DisplayName', 'High-end TEC Control (5 mK)');
xline(20, 'r-.', 'LineWidth', 1.4, 'DisplayName', 'Standard Lab TEC (20 mK)');
xline(50, 'g-.', 'LineWidth', 1.4, 'DisplayName', 'Open-loop Driver (50 mK)');

grid on; box on;
xlim([1, 100]);
ylim([min(nmse_lower_dB) - 2, max(nmse_upper_dB) + 2]);
xlabel('Thermal Control Accuracy \sigma_T [mK]');
ylabel('Normalized Mean Square Error (NMSE) [dB]');
title('Solver Accuracy Degradation vs. On-Chip Thermal Control Precision');
legend('Location', 'southeast');

%% 10. Export Figures to PNG
fprintf('Saving figures as PNG images...\n');
saveas(fig1, 'fig1_time_domain_thermal.png');
saveas(fig2, 'fig2_baseband_spectral_thermal.png');
saveas(fig3, 'fig3_monte_carlo_thermal.png');
saveas(fig4, 'fig4_sensitivity_sweep_thermal.png');
fprintf('Figures saved successfully in current directory.\n');
fprintf('========================================================================\n');
fprintf('Thermal error simulation complete!\n');
fprintf('========================================================================\n');
