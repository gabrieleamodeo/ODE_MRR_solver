

clear; close all; clc;

%% 1. Nominal Physical Parameters (Wu et al. 2014)
c             = 2.99792458e8;       % Speed of light in vacuum [m/s]
lambda0       = 1550.391e-9;        % Resonant wavelength [m] (1550.391 nm)
f0            = c / lambda0;        % Optical carrier frequency [Hz]
w0            = 2 * pi * f0;        % Optical carrier angular frequency [rad/s]
L             = 178.98e-6;          % Ring circumference [m] (178.98 um)
ng            = 4.1850;             % Waveguide group index
tau_rt        = ng * L / c;         % Cavity round-trip time [s] (~2.498 ps)
FSR           = 1 / tau_rt;         % Free Spectral Range [Hz] (~400.24 GHz)

% Target couplers (fixed microheaters)
kappa1        = 0.08;               % Coupler 1 (bus/through port power coupling)
kappa2        = 0.04;               % Coupler 2 (drop port power coupling)
r1            = sqrt(1 - kappa1);   % Field transmission coefficient coupler 1
r2            = sqrt(1 - kappa2);   % Field transmission coefficient coupler 2

% Nominal waveguide loss
alpha_nom_dB  = 8.0;                % Nominal waveguide power loss [dB/cm]
alpha_nom_lin = alpha_nom_dB * 100 * log(10) / 10; % [1/m]
eta_nom       = 1 - exp(-alpha_nom_lin * L);       % Nominal round-trip power loss
a_loss_nom    = sqrt(1 - eta_nom);                 % Nominal round-trip amplitude factor

%% 2. Coupled-Mode Theory (CMT) Formulation & Nominal Ideal Target ODE
% Cavity quality factors:
% Qi  = - (w0 * ng * L) / (c * ln(1 - eta))
% Qe1 = - (w0 * ng * L) / (c * ln(1 - kappa1))
% Qe2 = - (w0 * ng * L) / (c * ln(1 - kappa2))
Qi_nom        = - (w0 * ng * L) / (c * log(1 - eta_nom));
Qe1           = - (w0 * ng * L) / (c * log(1 - kappa1));
Qe2           = - (w0 * ng * L) / (c * log(1 - kappa2));

% Decay rates (gamma = w0 / (2*Q)):
gamma_i_nom   = w0 / (2 * Qi_nom);   % Nominal intrinsic decay rate [rad/s]
gamma_e1      = w0 / (2 * Qe1);      % Coupling decay rate to bus waveguide [rad/s]
gamma_e2      = w0 / (2 * Qe2);      % Coupling decay rate to drop waveguide [rad/s]

% Nominal ODE coefficients:
% Target ODE: dy(t)/dt + a0_nom*y(t) = dx(t)/dt + b0_nom*x(t)
a0_nom        = gamma_i_nom + gamma_e1 + gamma_e2;
b0_nom        = gamma_i_nom + gamma_e2 - gamma_e1;

%% 3. Analytical Critical Coupling Condition
% At critical coupling: b0 = 0 <=> gamma_i = gamma_e1 - gamma_e2
% Exactly corresponds to through-port field cancellation: r1 = r2 * a_loss
alpha_crit_lin   = - (1 / L) * log((1 - kappa1) / (1 - kappa2));
alpha_crit_dB_m  = alpha_crit_lin * 10 / log(10);
alpha_crit_dB_cm = alpha_crit_dB_m / 100;

fprintf('========================================================================\n');
fprintf('     Silicon Photonic ODE Solver - Cavity Propagation Loss Analysis     \n');
fprintf('========================================================================\n');
fprintf('Resonant Wavelength (lambda0)       : %.3f nm\n', lambda0 * 1e9);
fprintf('Optical Carrier Frequency (f0)      : %.3f THz\n', f0 / 1e12);
fprintf('Ring Circumference (L)              : %.2f um\n', L * 1e6);
fprintf('Waveguide Group Index (ng)          : %.4f\n', ng);
fprintf('Round-trip Time (tau_rt)            : %.3f ps\n', tau_rt * 1e12);
fprintf('Free Spectral Range (FSR)           : %.2f GHz\n', FSR / 1e9);
fprintf('Coupler 1 Power Coupling (kappa1)   : %.4f (r1 = %.4f)\n', kappa1, r1);
fprintf('Coupler 2 Power Coupling (kappa2)   : %.4f (r2 = %.4f)\n', kappa2, r2);
fprintf('Nominal Waveguide Loss (alpha_nom)  : %.2f dB/cm (eta = %.4f)\n', alpha_nom_dB, eta_nom);
fprintf('Critical Waveguide Loss (alpha_crit): %.4f dB/cm\n', alpha_crit_dB_cm);
fprintf('------------------------------------------------------------------------\n');
fprintf('Nominal Quality Factors             : Qi = %.2e, Qe1 = %.2e, Qe2 = %.2e\n', Qi_nom, Qe1, Qe2);
fprintf('Nominal Decay Rates                 : gamma_i = %.3e, gamma_e1 = %.3e, gamma_e2 = %.3e rad/s\n', ...
    gamma_i_nom, gamma_e1, gamma_e2);
fprintf('Nominal ODE Coefficients            : a0_nom = %.4e rad/s, b0_nom = %+.4e rad/s\n', a0_nom, b0_nom);
fprintf('Nominal Coupling Regime             : Over-coupled (b0_nom < 0)\n');
fprintf('========================================================================\n\n');

%% 4. Signal Definition: 10 Gb/s Gaussian Pulse (FWHM = 45 ps)
FWHM       = 45e-12;             % Full Width at Half Maximum: 45 ps
t0         = 0.3e-9;             % Pulse center: 300 ps
t_start    = -0.5e-9;            % -0.5 ns
t_end      = 1.5e-9;             % 1.5 ns
N          = 65536;              % Power of 2 grid points
t          = linspace(t_start, t_end, N);
dt         = t(2) - t(1);
fs         = 1 / dt;
f          = linspace(-fs/2, fs/2 - fs/N, N); % Baseband frequency axis [Hz]
w          = 2 * pi * f;                      % Baseband angular frequency [rad/s]

% Gaussian input pulse envelope
x_fun      = @(t_val) exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);
x          = x_fun(t);
I_in       = abs(x).^2;
E_in       = sum(I_in) * dt;                  % Total input pulse energy
X_w        = fftshift(fft(x));                % Input spectrum

% Ideal Target ODE Output (Nominal Reference)
% T_ODE(w) = (j*w + b0_nom) / (j*w + a0_nom)
T_ode_nom    = (1j * w + b0_nom) ./ (1j * w + a0_nom);
Y_ode_nom_w  = X_w .* T_ode_nom;
y_ode_nom    = ifft(ifftshift(Y_ode_nom_w));
I_ode_nom    = abs(y_ode_nom).^2;
E_ode_nom    = sum(I_ode_nom) * dt;

%% 5. Cavity Loss Parameter Sweep [4.0 to 14.0] dB/cm
N_sweep      = 101; % 101 points (0.1 dB/cm resolution)
alpha_vec    = linspace(4.0, 14.0, N_sweep);

% Pre-allocate metric arrays
a0_vec             = zeros(N_sweep, 1);
b0_vec             = zeros(N_sweep, 1);
notch_cmt_dB       = zeros(N_sweep, 1);
notch_mrr_dB       = zeros(N_sweep, 1);
NMSE_vec           = zeros(N_sweep, 1);
NMSE_dB_vec        = zeros(N_sweep, 1);
through_energy_pct = zeros(N_sweep, 1);
drop_leakage_pct   = zeros(N_sweep, 1);
drop_leakage_dB    = zeros(N_sweep, 1);
dissipated_pct     = zeros(N_sweep, 1);

fprintf('Sweeping waveguide loss across [%.1f, %.1f] dB/cm (%d points)...\n', ...
    alpha_vec(1), alpha_vec(end), N_sweep);

for i = 1:N_sweep
    a_dB_cm   = alpha_vec(i);
    a_dB_m    = a_dB_cm * 100;
    a_lin     = a_dB_m * log(10) / 10;
    eta_i     = 1 - exp(-a_lin * L);
    a_loss_i  = sqrt(1 - eta_i);
    
    % Quality factors and decay rates
    Qi_i      = - (w0 * ng * L) / (c * log(1 - eta_i));
    gamma_i   = w0 / (2 * Qi_i);
    
    % Realized ODE coefficients
    a0_i      = gamma_i + gamma_e1 + gamma_e2;
    b0_i      = gamma_i + gamma_e2 - gamma_e1;
    a0_vec(i) = a0_i;
    b0_vec(i) = b0_i;
    
    % Static Notch Depth / Extinction Ratio
    notch_cmt_dB(i) = 10 * log10((b0_i / a0_i)^2);
    notch_mrr_dB(i) = 10 * log10(abs((r1 - r2 * a_loss_i) / (1 - r1 * r2 * a_loss_i))^2);
    
    % Transfer functions
    % Through-port:
    H_through = (r1 - r2 * a_loss_i * exp(-1j * w * tau_rt)) ./ ...
                (1 - r1 * r2 * a_loss_i * exp(-1j * w * tau_rt));
    % Drop-port (optical leakage channel):
    H_drop    = (-sqrt(kappa1 * kappa2 * a_loss_i) * exp(-1j * w * tau_rt / 2)) ./ ...
                (1 - r1 * r2 * a_loss_i * exp(-1j * w * tau_rt));
    
    % Time-domain outputs via FFT/IFFT
    y_through = ifft(ifftshift(X_w .* H_through));
    y_drop    = ifft(ifftshift(X_w .* H_drop));
    
    I_through = abs(y_through).^2;
    I_drop    = abs(y_drop).^2;
    
    % Energy calculations
    E_through = sum(I_through) * dt;
    E_drop    = sum(I_drop) * dt;
    E_diss    = max(0, E_in - (E_through + E_drop));
    
    through_energy_pct(i) = (E_through / E_in) * 100;
    drop_leakage_pct(i)   = (E_drop / E_in) * 100;
    drop_leakage_dB(i)    = 10 * log10(E_drop / E_in);
    dissipated_pct(i)     = (E_diss / E_in) * 100;
    
    % NMSE w.r.t nominal ideal target ODE
    nmse_val       = sum((I_through - I_ode_nom).^2) / sum(I_ode_nom.^2);
    NMSE_vec(i)    = nmse_val;
    NMSE_dB_vec(i) = 10 * log10(nmse_val);
end

fprintf('Sweep complete.\n\n');

%% 6. Evaluation of 3 Representative Operating Points
% Point A: Under-loss (alpha = 5.0 dB/cm, deep over-coupling)
% Point B: Critical coupling (alpha = alpha_crit ~ 10.33 dB/cm, differentiator mode)
% Point C: Over-loss (alpha = 12.0 dB/cm, under-coupling)

alpha_pts = [5.0, alpha_crit_dB_cm, 12.0];
pt_names  = {'Point A: Under-Loss (Over-Coupled)', ...
             'Point B: Critical Coupling (Differentiator)', ...
             'Point C: Over-Loss (Under-Coupled)'};

y_through_pts = zeros(N, 3);
y_drop_pts    = zeros(N, 3);
I_through_pts = zeros(N, 3);
I_drop_pts    = zeros(N, 3);
metrics_pts   = struct();

fprintf('========================================================================\n');
fprintf('               Representative Operating Points Summary                  \n');
fprintf('========================================================================\n');

for p = 1:3
    a_p       = alpha_pts(p);
    a_p_lin   = a_p * 100 * log(10) / 10;
    eta_p     = 1 - exp(-a_p_lin * L);
    a_loss_p  = sqrt(1 - eta_p);
    
    Qi_p      = - (w0 * ng * L) / (c * log(1 - eta_p));
    gamma_i_p = w0 / (2 * Qi_p);
    a0_p      = gamma_i_p + gamma_e1 + gamma_e2;
    b0_p      = gamma_i_p + gamma_e2 - gamma_e1;
    
    H_thr_p   = (r1 - r2 * a_loss_p * exp(-1j * w * tau_rt)) ./ ...
                (1 - r1 * r2 * a_loss_p * exp(-1j * w * tau_rt));
    H_drp_p   = (-sqrt(kappa1 * kappa2 * a_loss_p) * exp(-1j * w * tau_rt / 2)) ./ ...
                (1 - r1 * r2 * a_loss_p * exp(-1j * w * tau_rt));
            
    y_thr_p   = ifft(ifftshift(X_w .* H_thr_p));
    y_drp_p   = ifft(ifftshift(X_w .* H_drp_p));
    
    y_through_pts(:, p) = y_thr_p;
    y_drop_pts(:, p)    = y_drp_p;
    I_through_pts(:, p) = abs(y_thr_p).^2;
    I_drop_pts(:, p)    = abs(y_drp_p).^2;
    
    E_thr_p   = sum(abs(y_thr_p).^2) * dt;
    E_drp_p   = sum(abs(y_drp_p).^2) * dt;
    E_diss_p  = max(0, E_in - (E_thr_p + E_drp_p));
    nmse_p    = sum((abs(y_thr_p).^2 - I_ode_nom).^2) / sum(I_ode_nom.^2);
    notch_p   = 10 * log10(abs((r1 - r2 * a_loss_p) / (1 - r1 * r2 * a_loss_p))^2);
    
    metrics_pts(p).alpha     = a_p;
    metrics_pts(p).a0        = a0_p;
    metrics_pts(p).b0        = b0_p;
    metrics_pts(p).notch_dB  = notch_p;
    metrics_pts(p).NMSE_dB   = 10 * log10(nmse_p);
    metrics_pts(p).E_thr_pct = (E_thr_p / E_in) * 100;
    metrics_pts(p).E_drp_pct = (E_drp_p / E_in) * 100;
    metrics_pts(p).E_dis_pct = (E_diss_p / E_in) * 100;
    
    fprintf('%s (alpha = %.2f dB/cm):\n', pt_names{p}, a_p);
    fprintf('  Coefficients        : a0 = %.4e rad/s, b0 = %+.4e rad/s\n', a0_p, b0_p);
    fprintf('  Resonance Notch     : %.2f dB\n', notch_p);
    fprintf('  NMSE vs Target ODE  : %.2f dB (Linear: %.4e)\n', 10*log10(nmse_p), nmse_p);
    fprintf('  Energy Budget       : Through = %.2f%%, Drop Leakage = %.2f%%, Cavity Loss = %.2f%%\n', ...
        metrics_pts(p).E_thr_pct, metrics_pts(p).E_drp_pct, metrics_pts(p).E_dis_pct);
    fprintf('------------------------------------------------------------------------\n');
end
fprintf('========================================================================\n\n');

%% 7. Visualization: Figure 1 - Metric Curves vs Waveguide Loss (alpha)
set(0, 'DefaultAxesFontSize', 10, 'DefaultLineLineWidth', 1.8);

fig1 = figure('Name', 'Cavity Loss Sweep - Solver Sensitivity & Metrics', ...
              'Position', [60, 60, 1050, 820], 'Color', 'w');

% (a) Realized ODE coefficients a0 and b0
subplot(2, 2, 1);
plot(alpha_vec, a0_vec / 1e10, 'b-', 'LineWidth', 2.0, 'DisplayName', 'ODE Coefficient a_0');
hold on;
plot(alpha_vec, b0_vec / 1e10, 'r-', 'LineWidth', 2.0, 'DisplayName', 'ODE Coefficient b_0');
yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Critical Coupling (b_0 = 0)');
xline(alpha_nom_dB, 'g-.', 'LineWidth', 1.4, 'DisplayName', sprintf('Nominal \\alpha (%.1f dB/cm)', alpha_nom_dB));
xline(alpha_crit_dB_cm, 'm:', 'LineWidth', 1.6, 'DisplayName', sprintf('Critical \\alpha_{crit} (%.2f dB/cm)', alpha_crit_dB_cm));
grid on; box on;
xlim([4.0, 14.0]);
xlabel('Waveguide Propagation Loss \alpha [dB/cm]');
ylabel('Coefficients [\times 10^{10} rad/s]');
title('(a) Realized ODE Coefficients vs. Waveguide Loss');
legend('Location', 'northwest');
% Annotate coupling regimes
text(5.0, 0.4, 'Over-coupled\n(b_0 < 0)', 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 9);
text(11.5, 0.4, 'Under-coupled\n(b_0 > 0)', 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 9);

% (b) Resonance Notch Depth / Extinction Ratio
subplot(2, 2, 2);
plot(alpha_vec, notch_mrr_dB, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Physical MRR Notch 10log_{10}(|H(0)|^2)');
hold on;
plot(alpha_vec, notch_cmt_dB, 'b--', 'LineWidth', 1.5, 'DisplayName', 'CMT Approx 10log_{10}(|b_0/a_0|^2)');
xline(alpha_crit_dB_cm, 'm:', 'LineWidth', 1.6, 'DisplayName', sprintf('\\alpha_{crit} = %.2f dB/cm', alpha_crit_dB_cm));
xline(alpha_nom_dB, 'g-.', 'LineWidth', 1.4, 'DisplayName', sprintf('Nominal \\alpha = %.1f dB/cm', alpha_nom_dB));
grid on; box on;
xlim([4.0, 14.0]);
ylim([-65, -10]);
xlabel('Waveguide Propagation Loss \alpha [dB/cm]');
ylabel('Resonance Notch Depth [dB]');
title('(b) Through-Port Static Extinction Ratio / Notch Depth');
legend('Location', 'south');

% (c) Through-Port NMSE relative to Nominal Ideal ODE Target
subplot(2, 2, 3);
plot(alpha_vec, NMSE_dB_vec, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Through-Port NMSE [dB]');
hold on;
xline(alpha_nom_dB, 'g-.', 'LineWidth', 1.4, 'DisplayName', sprintf('Nominal \\alpha (%.1f dB/cm)', alpha_nom_dB));
xline(alpha_crit_dB_cm, 'm:', 'LineWidth', 1.6, 'DisplayName', sprintf('\\alpha_{crit} (%.2f dB/cm)', alpha_crit_dB_cm));
plot(metrics_pts(1).alpha, metrics_pts(1).NMSE_dB, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Point A (5 dB/cm)');
plot(metrics_pts(2).alpha, metrics_pts(2).NMSE_dB, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Point B (10.33 dB/cm)');
plot(metrics_pts(3).alpha, metrics_pts(3).NMSE_dB, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Point C (12 dB/cm)');
grid on; box on;
xlim([4.0, 14.0]);
xlabel('Waveguide Propagation Loss \alpha [dB/cm]');
ylabel('NMSE [dB] vs. Target ODE');
title('(c) Through-Port Solver Error Degradation (NMSE)');
legend('Location', 'north');

% (d) Optical Pulse Energy Budget Breakdown
subplot(2, 2, 4);
plot(alpha_vec, through_energy_pct, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Through-Port Transmission (%)');
hold on;
plot(alpha_vec, drop_leakage_pct, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Drop-Port Leakage (%)');
plot(alpha_vec, dissipated_pct, 'k-.', 'LineWidth', 1.8, 'DisplayName', 'Cavity Dissipated / Absorbed (%)');
xline(alpha_nom_dB, 'g-.', 'LineWidth', 1.2, 'DisplayName', 'Nominal \alpha');
xline(alpha_crit_dB_cm, 'm:', 'LineWidth', 1.2, 'DisplayName', '\alpha_{crit}');
grid on; box on;
xlim([4.0, 14.0]);
ylim([0, 100]);
xlabel('Waveguide Propagation Loss \alpha [dB/cm]');
ylabel('Pulse Energy Fraction [%]');
title('(d) Optical Pulse Energy Budget Breakdown');
legend('Location', 'east');

%% 8. Visualization: Figure 2 - Time-Domain Waveforms at Operating Points A, B, C
fig2 = figure('Name', 'Cavity Loss Sweep - Time Domain Waveforms', ...
              'Position', [100, 100, 1000, 850], 'Color', 'w');

t_ps = t * 1e12;
t_zoom = [150, 650]; % Zoom around pulse center (300 ps)

for p = 1:3
    subplot(3, 1, p);
    plot(t_ps, I_in, 'k:', 'LineWidth', 1.4, 'DisplayName', 'Input Pulse |x(t)|^2');
    hold on;
    plot(t_ps, I_ode_nom, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Nominal Ideal ODE |y_{ODE}(t)|^2');
    plot(t_ps, I_through_pts(:, p), 'r-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Through-Port |y_{through}(t)|^2 (NMSE = %.1f dB)', metrics_pts(p).NMSE_dB));
    plot(t_ps, I_drop_pts(:, p), 'm-.', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('Drop-Port Leakage |y_{drop}(t)|^2 (Energy = %.1f%%)', metrics_pts(p).E_drp_pct));
    
    grid on; box on;
    xlim(t_zoom);
    xlabel('Time [ps]');
    ylabel('Intensity [a.u.]');
    
    switch p
        case 1
            title(sprintf('(a) Point A: Under-Loss / Deep Over-Coupled (\\alpha = %.1f dB/cm, b_0 = %+.2e rad/s)', ...
                metrics_pts(p).alpha, metrics_pts(p).b0));
        case 2
            title(sprintf('(b) Point B: Critical Coupling / Pure Differentiator (\\alpha = %.2f dB/cm, b_0 \\approx 0 rad/s)', ...
                metrics_pts(p).alpha));
        case 3
            title(sprintf('(c) Point C: Over-Loss / Under-Coupled (\\alpha = %.1f dB/cm, b_0 = %+.2e rad/s)', ...
                metrics_pts(p).alpha, metrics_pts(p).b0));
    end
    legend('Location', 'northeast');
end

%% 9. Export High-Resolution PNG Figures
fprintf('Saving figures as PNG images...\n');
saveas(fig1, 'fig1_loss_metrics_sweep.png');
saveas(fig2, 'fig2_loss_time_domain_waveforms.png');
fprintf('Figures saved successfully:\n');
fprintf(' - fig1_loss_metrics_sweep.png\n');
fprintf(' - fig2_loss_time_domain_waveforms.png\n');
fprintf('========================================================================\n');
fprintf('Cavity loss sensitivity analysis complete!\n');
fprintf('========================================================================\n');
