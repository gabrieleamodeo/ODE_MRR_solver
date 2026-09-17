
clear; close all; clc;

%% 1. Baseline Physical Constants & Waveguide Parameters (Wu et al. 2014)
c             = 2.99792458e8;       % Speed of light in vacuum [m/s]
lambda0       = 1550.391e-9;        % Target resonant carrier wavelength [m] (1550.391 nm)
f0            = c / lambda0;        % Optical carrier frequency [Hz] (~193.366 THz)
w0            = 2 * pi * f0;        % Optical carrier angular frequency [rad/s]
L             = 178.98e-6;          % Ring circumference [m] (178.98 um)
ng            = 4.1850;             % Silicon waveguide group index (TE mode)
tau_rt        = ng * L / c;         % Cavity round-trip time [s] (~2.498 ps)
FSR           = 1 / tau_rt;         % Free Spectral Range [Hz] (~400.24 GHz)

% Waveguide Propagation Loss from Paper
alpha_dB_cm   = 8.0;                % Waveguide power loss [dB/cm]
alpha_dB_m    = alpha_dB_cm * 100;  % [dB/m]
alpha_lin     = alpha_dB_m * log(10) / 10; % Linear power attenuation [1/m]
eta           = 1 - exp(-alpha_lin * L);   % Round-trip intrinsic power loss (~3.25%)
a_loss        = sqrt(1 - eta);             % Round-trip field amplitude transmission
gamma_i       = c * alpha_lin / (2 * ng);  % Cavity intrinsic decay rate [rad/s] (~6.598e9 rad/s)
Qi            = w0 / (2 * gamma_i);        % Intrinsic cavity quality factor

% MZI Tunable Directional Coupler Parameters
kappa0        = 0.0441;             % Directional coupler splitting power ratio (FDTD)
kappa_max     = 4 * kappa0 * (1 - kappa0); % Max MZI power coupling (~0.168621)
dn_dT         = 1.86e-4;            % Thermo-optic coefficient of silicon [1/K]

% Physical lengths of bus arm microheaters (Wu et al. 2014):
% For MRR 1: Heater 11 (Lb11 = 116.80 um), Heater 12 (Lb12 = 47.12 um)
% For MRR 2: Heater 21 (Lb21 = 116.80 um), Heater 22 (Lb22 = 47.12 um)
Lb1           = 116.80e-6;          % Bus arm 1 heater length [m] (Through coupler)
Lb2           = 47.12e-6;           % Bus arm 2 heater length [m] (Drop coupler)

%% 2. Target Design Coefficients for the Cascaded Second-Order ODE
% Target power coupling coefficients for MRR 1 and MRR 2:
% To realize a general second-order ODE with distinct poles and non-zero zeros
% (and prevent the degenerate critically-damped pure-differentiator artifact
% where the carrier notch is quadratically flattened), we synthesize distinct,
% physically realizable coupling rates:
%   MRR 1 (Over-coupled stage) : kappa11 = 0.12, kappa12 = 0.02
%   MRR 2 (Under-coupled stage): kappa21 = 0.04, kappa22 = 0.06
kappa11_tgt   = 0.12;               % Through-port coupler power coupling MRR 1
kappa12_tgt   = 0.02;               % Drop-port coupler power coupling MRR 1
kappa21_tgt   = 0.04;               % Through-port coupler power coupling MRR 2
kappa22_tgt   = 0.06;               % Drop-port coupler power coupling MRR 2

% Nominal MZI phase differences needed to synthesize target couplers:
cos_dphi11    = (kappa11_tgt / (2 * kappa0 * (1 - kappa0))) - 1;
cos_dphi12    = (kappa12_tgt / (2 * kappa0 * (1 - kappa0))) - 1;
cos_dphi21    = (kappa21_tgt / (2 * kappa0 * (1 - kappa0))) - 1;
cos_dphi22    = (kappa22_tgt / (2 * kappa0 * (1 - kappa0))) - 1;

dphi11_nom    = acos(min(max(cos_dphi11, -1), 1)); % Nominal phase difference Coupler 11 [rad]
dphi12_nom    = acos(min(max(cos_dphi12, -1), 1)); % Nominal phase difference Coupler 12 [rad]
dphi21_nom    = acos(min(max(cos_dphi21, -1), 1)); % Nominal phase difference Coupler 21 [rad]
dphi22_nom    = acos(min(max(cos_dphi22, -1), 1)); % Nominal phase difference Coupler 22 [rad]

% Quality factors and decay rates of individual rings:
% Stage 1 (MRR 1):
Qe11_nom      = - (w0 * ng * L) / (c * log(1 - kappa11_tgt));
Qe12_nom      = - (w0 * ng * L) / (c * log(1 - kappa12_tgt));
gamma_e11_nom = w0 / (2 * Qe11_nom); % ~2.59e10 rad/s
gamma_e12_nom = w0 / (2 * Qe12_nom); % ~4.04e9 rad/s
a10_target    = gamma_i + gamma_e11_nom + gamma_e12_nom; % ~3.62e10 rad/s
b10_target    = gamma_i + gamma_e12_nom - gamma_e11_nom; % ~-1.49e10 rad/s

% Stage 2 (MRR 2: Distinct pole and zero):
Qe21_nom      = - (w0 * ng * L) / (c * log(1 - kappa21_tgt));
Qe22_nom      = - (w0 * ng * L) / (c * log(1 - kappa22_tgt));
gamma_e21_nom = w0 / (2 * Qe21_nom); % ~8.17e9 rad/s
gamma_e22_nom = w0 / (2 * Qe22_nom); % ~1.24e10 rad/s
a20_target    = gamma_i + gamma_e21_nom + gamma_e22_nom; % ~2.72e10 rad/s
b20_target    = gamma_i + gamma_e22_nom - gamma_e21_nom; % ~+1.08e10 rad/s

% Cascaded Second-Order ODE Coefficients (Wu et al. Eq. 9):
%   d^2y/dt^2 + a1 * dy/dt + a0 * y = b2 * d^2x/dt^2 + b1 * dx/dt + b0 * x
a1_target     = a10_target + a20_target;               % [rad/s] (~6.34e10 rad/s)
a0_target     = a10_target * a20_target;               % [(rad/s)^2] (~9.83e20 (rad/s)^2)
b2_target     = 1.0;                                   % Dimensionless
b1_target     = b10_target + b20_target;               % [rad/s] (~-4.13e9 rad/s)
b0_target     = b10_target * b20_target;               % [(rad/s)^2] (~-1.62e20 (rad/s)^2)

%% 3. Thermal Noise & Control Accuracy Settings
% Standard laboratory temperature control precision (Standard Deviation sigma_T):
% You can adjust sigma_T to match your lab equipment!
sigma_T_heater= 0.020;              % [K] 20 mK standard deviation of heater drivers
sigma_T_ring  = 0.015;              % [K] 15 mK residual cavity substrate drift
thermal_crosstalk = 0.06;           % 6% thermal leakage from microheaters to ring cavity

fprintf('========================================================================================\n');
fprintf('     SECOND-ORDER OPTICAL ODE SOLVER VIA CASCADED SILICON PHOTONIC MRRs                 \n');
fprintf('     Reference: J. Wu et al., Optics Express 22(21), 26254-26264 (2014) [Eq. (9)]       \n');
fprintf('========================================================================================\n');
fprintf('Waveguide Propagation Loss (alpha)  : %.2f dB/cm (Intrinsic Qi = %.2e)\n', alpha_dB_cm, Qi);
fprintf('Stage 1 Target Coefficients         : a10 = %.4e rad/s, b10 = %+.4e rad/s\n', a10_target, b10_target);
fprintf('Stage 2 Target Coefficients         : a20 = %.4e rad/s, b20 = %+.4e rad/s\n', a20_target, b20_target);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('SYNTHESIZED SECOND-ORDER ODE (d^2y/dt^2 + a1*dy/dt + a0*y = d^2x/dt^2 + b1*dx/dt + b0*x):\n');
fprintf('  Damping / Bandwidth Rate (a1)     : %.4e rad/s\n', a1_target);
fprintf('  Resonance Stiffness Term (a0)     : %.4e (rad/s)^2\n', a0_target);
fprintf('  First-Derivative Zero Weight (b1) : %+.4e rad/s\n', b1_target);
fprintf('  Zero Constant Coefficient (b0)    : %+.4e (rad/s)^2\n', b0_target);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('INDEPENDENT THERMAL CONTROL INACCURACIES:\n');
fprintf('  Heater Driver Jitter (sigma_T)    : %.2f mK (4 Independent Heaters: H11, H12, H21, H22)\n', sigma_T_heater * 1e3);
fprintf('  Substrate Drift / Cavity Jitter   : %.2f mK (2 Independent Cavities: Ring 1, Ring 2)\n', sigma_T_ring * 1e3);
fprintf('========================================================================================\n\n');

%% 4. Signal Definition: 10 Gb/s Gaussian Pulse (FWHM = 45 ps)
FWHM          = 45e-12;             % 45 ps
t0            = 0.3e-9;             % Center: 300 ps
t_start       = -0.5e-9;
t_end         = 1.5e-9;
N             = 65536;              % Power of 2 grid points
t             = linspace(t_start, t_end, N);
dt            = t(2) - t(1);
fs            = 1 / dt;
f             = linspace(-fs/2, fs/2 - fs/N, N); % Baseband frequency axis [Hz]
w             = 2 * pi * f;                      % Baseband angular frequency [rad/s]

% Input Gaussian pulse envelope
x_fun         = @(t_val) exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);
x             = x_fun(t);
I_in          = abs(x).^2;
X_w           = fftshift(fft(x));

% -------------------------------------------------------------------------
% (A) Ideal Target Transfer Functions (Coupled-Mode Theory):
% -------------------------------------------------------------------------
% 1st-Order Target ODE Transfer Function:
T_ode_1st     = (1j * w + b10_target) ./ (1j * w + a10_target);
y_ode_1st     = ifft(ifftshift(X_w .* T_ode_1st));
I_ode_1st     = abs(y_ode_1st).^2;

% 2nd-Order Target ODE Transfer Function (Eq. 9 in paper):
% T_ODE,2(w) = [(jw)^2 + b1*(jw) + b0] / [(jw)^2 + a1*(jw) + a0]
T_ode_2nd     = ((1j * w).^2 + b1_target * (1j * w) + b0_target) ./ ...
                ((1j * w).^2 + a1_target * (1j * w) + a0_target);
y_ode_2nd     = ifft(ifftshift(X_w .* T_ode_2nd));
I_ode_2nd     = abs(y_ode_2nd).^2;

% -------------------------------------------------------------------------
% (B) Nominal Physical Cascaded MRR Model (Noiseless):
% -------------------------------------------------------------------------
r11_nom       = sqrt(1 - kappa11_tgt);
r12_nom       = sqrt(1 - kappa12_tgt);
r21_nom       = sqrt(1 - kappa21_tgt);
r22_nom       = sqrt(1 - kappa22_tgt);

H1_nom        = (r11_nom - r12_nom * a_loss * exp(-1j * w * tau_rt)) ./ ...
                (1 - r11_nom * r12_nom * a_loss * exp(-1j * w * tau_rt));
H2_nom        = (r21_nom - r22_nom * a_loss * exp(-1j * w * tau_rt)) ./ ...
                (1 - r21_nom * r22_nom * a_loss * exp(-1j * w * tau_rt));
H_casc_nom    = H1_nom .* H2_nom;

y_casc_nom    = ifft(ifftshift(X_w .* H_casc_nom));
I_casc_nom    = abs(y_casc_nom).^2;

NMSE_casc_nom = sum((I_ode_2nd - I_casc_nom).^2) / sum(I_ode_2nd.^2);
NMSE_1st_nom  = sum((I_ode_1st - abs(ifft(ifftshift(X_w .* H1_nom))).^2).^2) / sum(I_ode_1st.^2);

%% 5. Single Perturbed Realization (Realistic Lab Shot with 4 Independent Heaters)
rng(42); % Fixed seed for reproducible single-shot display

% Sample 4 independent heater temperature errors:
dT_h11        = randn * sigma_T_heater;
dT_h12        = randn * sigma_T_heater;
dT_h21        = randn * sigma_T_heater;
dT_h22        = randn * sigma_T_heater;

% Sample 2 independent cavity substrate drifts + heater cross-talk:
dT_cav1       = randn * sigma_T_ring + thermal_crosstalk * (dT_h11 + dT_h12) / 2;
dT_cav2       = randn * sigma_T_ring + thermal_crosstalk * (dT_h21 + dT_h22) / 2;

% Thermo-optic phase perturbations for Ring 1:
dphi11_noisy  = dphi11_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dT_h11;
dphi12_noisy  = dphi12_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dT_h12;
kappa11_noisy = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dphi11_noisy)), 1e-6), 0.99);
kappa12_noisy = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dphi12_noisy)), 1e-6), 0.99);
r11_n         = sqrt(1 - kappa11_noisy);
r12_n         = sqrt(1 - kappa12_noisy);
dphi_rt1_n    = (2 * pi / lambda0) * dn_dT * L * dT_cav1;
Delta_f_res1  = dphi_rt1_n / (2 * pi * tau_rt);

% Thermo-optic phase perturbations for Ring 2 (Independent!):
dphi21_noisy  = dphi21_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dT_h21;
dphi22_noisy  = dphi22_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dT_h22;
kappa21_noisy = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dphi21_noisy)), 1e-6), 0.99);
kappa22_noisy = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dphi22_noisy)), 1e-6), 0.99);
r21_n         = sqrt(1 - kappa21_noisy);
r22_n         = sqrt(1 - kappa22_noisy);
dphi_rt2_n    = (2 * pi / lambda0) * dn_dT * L * dT_cav2;
Delta_f_res2  = dphi_rt2_n / (2 * pi * tau_rt);

% Inter-cavity detuning between Ring 1 and Ring 2:
Delta_f_inter = abs(Delta_f_res1 - Delta_f_res2);

% Transfer functions of the two noisy rings:
H1_noisy      = (r11_n - r12_n * a_loss * exp(-1j * (w * tau_rt + dphi_rt1_n))) ./ ...
                (1 - r11_n * r12_n * a_loss * exp(-1j * (w * tau_rt + dphi_rt1_n)));
H2_noisy      = (r21_n - r22_n * a_loss * exp(-1j * (w * tau_rt + dphi_rt2_n))) ./ ...
                (1 - r21_n * r22_n * a_loss * exp(-1j * (w * tau_rt + dphi_rt2_n)));
H_casc_noisy  = H1_noisy .* H2_noisy;

% Time-domain output intensity:
y_casc_noisy  = ifft(ifftshift(X_w .* H_casc_noisy));
I_casc_noisy  = abs(y_casc_noisy).^2;
NMSE_casc_noisy = sum((I_ode_2nd - I_casc_noisy).^2) / sum(I_ode_2nd.^2);

% Realized 2nd-order ODE coefficients:
ge11_n        = - c * log(1 - kappa11_noisy) / (2 * ng * L);
ge12_n        = - c * log(1 - kappa12_noisy) / (2 * ng * L);
ge21_n        = - c * log(1 - kappa21_noisy) / (2 * ng * L);
ge22_n        = - c * log(1 - kappa22_noisy) / (2 * ng * L);

a10_n         = gamma_i + ge11_n + ge12_n;
b10_n         = gamma_i + ge12_n - ge11_n;
a20_n         = gamma_i + ge21_n + ge22_n;
b20_n         = gamma_i + ge22_n - ge21_n;

a1_realized   = a10_n + a20_n;
a0_realized   = a10_n * a20_n;
b1_realized   = b10_n + b20_n;
b0_realized   = b10_n * b20_n;

fprintf('--- SINGLE REALIZATION BREAKDOWN (Independent Heaters) ---\n');
fprintf('MRR 1: dT11 = %+.2f mK (k11=%.4f), dT12 = %+.2f mK (k12=%.4f)\n', ...
        dT_h11*1e3, kappa11_noisy, dT_h12*1e3, kappa12_noisy);
fprintf('       Resonance Detuning 1: %+.1f MHz\n', Delta_f_res1 / 1e6);
fprintf('MRR 2: dT21 = %+.2f mK (k21=%.4f), dT22 = %+.2f mK (k22=%.4f)\n', ...
        dT_h21*1e3, kappa21_noisy, dT_h22*1e3, kappa22_noisy);
fprintf('       Resonance Detuning 2: %+.1f MHz\n', Delta_f_res2 / 1e6);
fprintf('Inter-Cavity Detuning (|df1 - df2|) : %.1f MHz (Causes notch splitting!)\n', Delta_f_inter / 1e6);
fprintf('Target ODE Coefficients             : a1 = %.4e rad/s, a0 = %.4e (rad/s)^2\n', a1_target, a0_target);
fprintf('                                      b1 = %+.4e rad/s, b0 = %+.4e (rad/s)^2\n', b1_target, b0_target);
fprintf('Realized Noisy Coefficients         : a1 = %.4e rad/s, a0 = %.4e (rad/s)^2\n', a1_realized, a0_realized);
fprintf('                                      b1 = %+.4e rad/s, b0 = %+.4e (rad/s)^2\n', b1_realized, b0_realized);
fprintf('Nominal Cascaded MRR NMSE (Noiseless): %.4e (%.2f dB)\n', NMSE_casc_nom, 10*log10(NMSE_casc_nom));
fprintf('Thermally Perturbed Cascaded NMSE   : %.4e (%.2f dB)\n', NMSE_casc_noisy, 10*log10(NMSE_casc_noisy));
fprintf('NMSE Degradation due to Noise       : %+.2f dB\n', 10*log10(NMSE_casc_noisy) - 10*log10(NMSE_casc_nom));
fprintf('========================================================================================\n\n');

%% 6. Monte Carlo Statistical Simulation (N_trials = 200 Realizations)
N_trials      = 200;
NMSE_mc_2nd   = zeros(N_trials, 1);
NMSE_mc_1st   = zeros(N_trials, 1);
I_mc_2nd      = zeros(N, N_trials);

a1_mc         = zeros(N_trials, 1);
a0_mc         = zeros(N_trials, 1);
b1_mc         = zeros(N_trials, 1);
b0_mc         = zeros(N_trials, 1);
Delta_f_inter_mc = zeros(N_trials, 1);

fprintf('Running Monte Carlo simulation (%d trials) with 4 independent heaters...\n', N_trials);
for trial = 1:N_trials
    % 4 independent heater fluctuations
    dt11 = randn * sigma_T_heater;
    dt12 = randn * sigma_T_heater;
    dt21 = randn * sigma_T_heater;
    dt22 = randn * sigma_T_heater;
    
    % 2 independent cavity fluctuations
    dt_c1 = randn * sigma_T_ring + thermal_crosstalk * (dt11 + dt12) / 2;
    dt_c2 = randn * sigma_T_ring + thermal_crosstalk * (dt21 + dt22) / 2;
    
    % Couplers MRR 1
    dp11 = dphi11_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt11;
    dp12 = dphi12_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt12;
    k11 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp11)), 1e-6), 0.99);
    k12 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp12)), 1e-6), 0.99);
    r11_t = sqrt(1 - k11); r12_t = sqrt(1 - k12);
    dp_rt1 = (2 * pi / lambda0) * dn_dT * L * dt_c1;
    
    % Couplers MRR 2
    dp21 = dphi21_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt21;
    dp22 = dphi22_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt22;
    k21 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp21)), 1e-6), 0.99);
    k22 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp22)), 1e-6), 0.99);
    r21_t = sqrt(1 - k21); r22_t = sqrt(1 - k22);
    dp_rt2 = (2 * pi / lambda0) * dn_dT * L * dt_c2;
    
    Delta_f_inter_mc(trial) = abs(dp_rt1 - dp_rt2) / (2 * pi * tau_rt);
    
    % Transfer functions
    H1_t = (r11_t - r12_t * a_loss * exp(-1j * (w * tau_rt + dp_rt1))) ./ ...
           (1 - r11_t * r12_t * a_loss * exp(-1j * (w * tau_rt + dp_rt1)));
    H2_t = (r21_t - r22_t * a_loss * exp(-1j * (w * tau_rt + dp_rt2))) ./ ...
           (1 - r21_t * r22_t * a_loss * exp(-1j * (w * tau_rt + dp_rt2)));
    H_casc_t = H1_t .* H2_t;
    
    y_casc_t = ifft(ifftshift(X_w .* H_casc_t));
    I_t = abs(y_casc_t).^2;
    I_mc_2nd(:, trial) = I_t;
    
    % NMSE for 2nd-order solver
    NMSE_mc_2nd(trial) = sum((I_ode_2nd - I_t).^2) / sum(I_ode_2nd.^2);
    
    % NMSE for 1st-order solver (Ring 1 alone for direct comparison)
    y1_t = ifft(ifftshift(X_w .* H1_t));
    NMSE_mc_1st(trial) = sum((I_ode_1st - abs(y1_t).^2).^2) / sum(I_ode_1st.^2);
    
    % Record realized coefficients
    ge11_t = - c * log(1 - k11) / (2 * ng * L);
    ge12_t = - c * log(1 - k12) / (2 * ng * L);
    ge21_t = - c * log(1 - k21) / (2 * ng * L);
    ge22_t = - c * log(1 - k22) / (2 * ng * L);
    a10_t  = gamma_i + ge11_t + ge12_t;
    b10_t  = gamma_i + ge12_t - ge11_t;
    a20_t  = gamma_i + ge21_t + ge22_t;
    b20_t  = gamma_i + ge22_t - ge21_t;
    a1_mc(trial) = a10_t + a20_t;
    a0_mc(trial) = a10_t * a20_t;
    b1_mc(trial) = b10_t + b20_t;
    b0_mc(trial) = b10_t * b20_t;
end

mean_NMSE_2nd = mean(NMSE_mc_2nd);
mean_NMSE_1st = mean(NMSE_mc_1st);
mean_I_2nd    = mean(I_mc_2nd, 2);
std_I_2nd     = std(I_mc_2nd, 0, 2);

fprintf('Monte Carlo Complete (%d Trials):\n', N_trials);
fprintf('  1st-Order Mean NMSE         : %.4e (%.2f dB)\n', mean_NMSE_1st, 10*log10(mean_NMSE_1st));
fprintf('  2nd-Order Cascaded Mean NMSE: %.4e (%.2f dB)\n', mean_NMSE_2nd, 10*log10(mean_NMSE_2nd));
fprintf('  2nd-Order Worst NMSE        : %.4e (%.2f dB)\n', max(NMSE_mc_2nd), 10*log10(max(NMSE_mc_2nd)));
fprintf('  2nd-Order Best NMSE         : %.4e (%.2f dB)\n', min(NMSE_mc_2nd), 10*log10(min(NMSE_mc_2nd)));
fprintf('========================================================================================\n\n');

%% 7. Thermal Sensitivity Sweep: NMSE vs. sigma_T (1st-Order vs. 2nd-Order)
fprintf('Running Sensitivity Sweep (NMSE vs. sigma_T across 11 precision tiers)...\n');
sigma_T_vec   = [1, 2, 5, 10, 15, 20, 30, 40, 50, 75, 100] * 1e-3; % 1 mK to 100 mK
N_sweep_pts   = length(sigma_T_vec);
sweep_trials  = 40;

NMSE_sweep_1st_dB = zeros(N_sweep_pts, 1);
NMSE_sweep_2nd_dB = zeros(N_sweep_pts, 1);

for s_idx = 1:N_sweep_pts
    sig_T = sigma_T_vec(s_idx);
    e1_arr = zeros(sweep_trials, 1);
    e2_arr = zeros(sweep_trials, 1);
    
    for tr = 1:sweep_trials
        dt11 = randn * sig_T; dt12 = randn * sig_T;
        dt21 = randn * sig_T; dt22 = randn * sig_T;
        dt_c1 = randn * (sig_T * 0.75) + thermal_crosstalk * (dt11 + dt12) / 2;
        dt_c2 = randn * (sig_T * 0.75) + thermal_crosstalk * (dt21 + dt22) / 2;
        
        dp11 = dphi11_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt11;
        dp12 = dphi12_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt12;
        k11 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp11)), 1e-6), 0.99);
        k12 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp12)), 1e-6), 0.99);
        
        dp21 = dphi21_nom + (2 * pi / lambda0) * dn_dT * Lb1 * dt21;
        dp22 = dphi22_nom + (2 * pi / lambda0) * dn_dT * Lb2 * dt22;
        k21 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp21)), 1e-6), 0.99);
        k22 = min(max(2 * kappa0 * (1 - kappa0) * (1 + cos(dp22)), 1e-6), 0.99);
        
        dp_rt1 = (2 * pi / lambda0) * dn_dT * L * dt_c1;
        dp_rt2 = (2 * pi / lambda0) * dn_dT * L * dt_c2;
        
        H1_s = (sqrt(1-k11) - sqrt(1-k12)*a_loss*exp(-1j*(w*tau_rt + dp_rt1))) ./ ...
               (1 - sqrt(1-k11)*sqrt(1-k12)*a_loss*exp(-1j*(w*tau_rt + dp_rt1)));
        H2_s = (sqrt(1-k21) - sqrt(1-k22)*a_loss*exp(-1j*(w*tau_rt + dp_rt2))) ./ ...
               (1 - sqrt(1-k21)*sqrt(1-k22)*a_loss*exp(-1j*(w*tau_rt + dp_rt2)));
        
        y1_s = ifft(ifftshift(X_w .* H1_s));
        y2_s = ifft(ifftshift(X_w .* (H1_s .* H2_s)));
        
        e1_arr(tr) = sum((I_ode_1st - abs(y1_s).^2).^2) / sum(I_ode_1st.^2);
        e2_arr(tr) = sum((I_ode_2nd - abs(y2_s).^2).^2) / sum(I_ode_2nd.^2);
    end
    NMSE_sweep_1st_dB(s_idx) = 10 * log10(mean(e1_arr));
    NMSE_sweep_2nd_dB(s_idx) = 10 * log10(mean(e2_arr));
end

fprintf('Sensitivity sweep complete.\n\n');

%% 8. Visualization: Figure 1 - Time-Domain Waveforms & Cascaded Behavior
set(0, 'DefaultAxesFontSize', 10.5, 'DefaultLineLineWidth', 1.8);

fig1 = figure('Name', '2nd-Order Cascaded ODE Solver - Time Domain', ...
              'Position', [50, 40, 1000, 850], 'Color', 'w');
t_ps = t * 1e12;

% (a) Input pulse vs 1st-order output vs 2nd-order target output
subplot(3, 1, 1);
plot(t_ps, I_in, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Input Gaussian Pulse |x(t)|^2');
hold on;
plot(t_ps, I_ode_1st, 'm--', 'LineWidth', 1.8, 'DisplayName', '1st-Order Target ODE (Single MRR)');
plot(t_ps, I_ode_2nd, 'b-', 'LineWidth', 2.2, 'DisplayName', '2nd-Order Target ODE (Cascaded MRRs)');
grid on; box on;
xlim([t_start*1e12, 1000]);
xlabel('Time [ps]'); ylabel('Normalized Intensity');
title('(a) Temporal Output Progression: Input Pulse \rightarrow 1st-Order \rightarrow 2nd-Order ODE');
legend('Location', 'northeast', 'FontSize', 9.0);

% (b) Ideal 2nd-order ODE vs Nominal MRR vs Thermally Perturbed MRR
subplot(3, 1, 2);
plot(t_ps, I_ode_2nd, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Ideal 2nd-Order ODE Target');
hold on;
plot(t_ps, I_casc_nom, 'g--', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Nominal Cascaded MRR (NMSE = %.1f dB)', 10*log10(NMSE_casc_nom)));
plot(t_ps, I_casc_noisy, 'r-.', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Cascaded MRR with Thermal Noise (\\sigma_T=20 mK, NMSE = %.1f dB)', 10*log10(NMSE_casc_noisy)));
grid on; box on;
xlim([t_start*1e12, 1000]);
xlabel('Time [ps]'); ylabel('Intensity [a.u.]');
title('(b) Second-Order Time-Domain Output: Ideal vs. Nominal vs. Thermally Perturbed MRRs');
legend('Location', 'northeast', 'FontSize', 9.0);

% (c) Instantaneous intensity errors
subplot(3, 1, 3);
plot(t_ps, abs(I_ode_2nd - I_casc_nom), 'g--', 'LineWidth', 1.5, ...
     'DisplayName', 'Error: |I_{ODE,2} - I_{Casc,nom}| (Intrinsic Periodicity Limit)');
hold on;
plot(t_ps, abs(I_ode_2nd - I_casc_noisy), 'r-', 'LineWidth', 1.8, ...
     'DisplayName', 'Error: |I_{ODE,2} - I_{Casc,noisy}| (4 Independent Heaters Jitter)');
grid on; box on;
xlim([t_start*1e12, 1000]);
xlabel('Time [ps]'); ylabel('Intensity Deviation');
title('(c) Instantaneous Deviation in 2nd-Order Output Due to Independent Heater Controls');
legend('Location', 'northeast', 'FontSize', 9.0);

%% 9. Visualization: Figure 2 - Baseband Spectral Response & Inter-Cavity Detuning
fig2 = figure('Name', '2nd-Order Cascaded ODE Solver - Spectral Response', ...
              'Position', [100, 80, 980, 680], 'Color', 'w');
f_GHz = f / 1e9;
zoom_GHz = 20; % +/- 20 GHz baseband zoom

% (a) Baseband Magnitude Response [dB]
subplot(2, 1, 1);
plot(f_GHz, 10 * log10(abs(T_ode_2nd).^2), 'b--', 'LineWidth', 2.2, 'DisplayName', 'Ideal 2nd-Order ODE Target');
hold on;
plot(f_GHz, 10 * log10(abs(H_casc_nom).^2), 'g-', 'LineWidth', 1.6, 'DisplayName', 'Nominal Cascaded MRR (Noiseless, Co-Resonant)');
plot(f_GHz, 10 * log10(abs(H_casc_noisy).^2), 'r-', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Perturbed Cascaded MRR (Inter-Cavity Detuning |\\Delta f| = %.1f MHz)', Delta_f_inter/1e6));
grid on; box on;
xlim([-zoom_GHz, zoom_GHz]); ylim([-50, 5]);
xlabel('Baseband Frequency Offset f - f_0 [GHz]');
ylabel('Transmission [dB]');
title('(a) Second-Order Resonance Notch: Deep Dual-Resonator Transmission & Detuning Splitting');
legend('Location', 'southeast', 'FontSize', 9.0);

% (b) Baseband Phase Response [deg]
subplot(2, 1, 2);
plot(f_GHz, angle(T_ode_2nd) * 180 / pi, 'b--', 'LineWidth', 2.2, 'DisplayName', 'Ideal 2nd-Order ODE Phase');
hold on;
plot(f_GHz, angle(H_casc_nom) * 180 / pi, 'g-', 'LineWidth', 1.6, 'DisplayName', 'Nominal Cascaded MRR');
plot(f_GHz, angle(H_casc_noisy) * 180 / pi, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Perturbed Cascaded MRR');
grid on; box on;
xlim([-zoom_GHz, zoom_GHz]); ylim([-190, 190]);
xlabel('Baseband Frequency Offset f - f_0 [GHz]');
ylabel('Phase [deg]');
title('(b) Second-Order Phase Response (Accumulated 360^\circ Phase Transition)');
legend('Location', 'southeast', 'FontSize', 9.0);

%% 10. Visualization: Figure 3 - Monte Carlo Statistical Analysis (N = 200 Trials)
fig3 = figure('Name', '2nd-Order Cascaded ODE Solver - Monte Carlo Statistics', ...
              'Position', [140, 120, 1050, 680], 'Color', 'w');

% (a) 2nd-Order Output Waveform Confidence Bands
subplot(2, 2, [1, 3]);
y_upper2 = (mean_I_2nd + 2*std_I_2nd)';
y_lower2 = (max(mean_I_2nd - 2*std_I_2nd, 0))';
y_upper1 = (mean_I_2nd + std_I_2nd)';
y_lower1 = (max(mean_I_2nd - std_I_2nd, 0))';

fill([t_ps, fliplr(t_ps)], [y_upper2, fliplr(y_lower2)], ...
     [1.0, 0.88, 0.88], 'EdgeColor', 'none', 'DisplayName', '\pm 2\sigma Band (95.4% Confidence)');
hold on;
fill([t_ps, fliplr(t_ps)], [y_upper1, fliplr(y_lower1)], ...
     [1.0, 0.70, 0.70], 'EdgeColor', 'none', 'DisplayName', '\pm 1\sigma Band (68.3% Confidence)');
plot(t_ps, I_ode_2nd, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Ideal 2nd-Order Target');
plot(t_ps, mean_I_2nd, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Monte Carlo Mean Waveform');
grid on; box on;
xlim([100, 600]);
xlabel('Time [ps]'); ylabel('Intensity [a.u.]');
title(sprintf('(a) Second-Order Waveform Uncertainty Bands (%d Trials, \\sigma_T = 20 mK)', N_trials));
legend('Location', 'northeast', 'FontSize', 9.0);

% (b) NMSE Distribution Comparison: 1st-Order vs 2nd-Order Cascaded
subplot(2, 2, 2);
histogram(10*log10(NMSE_mc_1st), 15, 'FaceColor', [0.3 0.6 0.9], 'FaceAlpha', 0.6, ...
          'DisplayName', sprintf('1st-Order Single MRR (Mean = %.1f dB)', 10*log10(mean_NMSE_1st)));
hold on;
histogram(10*log10(NMSE_mc_2nd), 15, 'FaceColor', [0.9 0.3 0.3], 'FaceAlpha', 0.6, ...
          'DisplayName', sprintf('2nd-Order Cascaded (Mean = %.1f dB)', 10*log10(mean_NMSE_2nd)));
xline(10*log10(mean_NMSE_1st), 'b--', 'LineWidth', 1.8, 'DisplayName', '1st-Order Mean');
xline(10*log10(mean_NMSE_2nd), 'r--', 'LineWidth', 1.8, 'DisplayName', '2nd-Order Mean');
grid on; box on;
xlabel('Through-Port NMSE [dB]'); ylabel('Occurrences');
title('(b) NMSE Statistical Distribution (1st-Order vs. 2nd-Order)');
legend('Location', 'northwest', 'FontSize', 8.5);

% (c) Inter-Cavity Detuning Histogram
subplot(2, 2, 4);
histogram(Delta_f_inter_mc / 1e6, 18, 'FaceColor', [0.5 0.2 0.7], 'FaceAlpha', 0.7, ...
          'DisplayName', 'Inter-Cavity Resonance Mismatch |\Delta f_{res,1} - \Delta f_{res,2}|');
xline(mean(Delta_f_inter_mc)/1e6, 'k--', 'LineWidth', 1.8, ...
      'DisplayName', sprintf('Mean Mismatch = %.1f MHz', mean(Delta_f_inter_mc)/1e6));
grid on; box on;
xlabel('Inter-Cavity Detuning |\Delta f_{res}| [MHz]'); ylabel('Occurrences');
title('(c) Inter-Cavity Resonance Jitter Between Cascaded MRRs');
legend('Location', 'northeast', 'FontSize', 8.5);

%% 11. Visualization: Figure 4 - Thermal Sensitivity Sweep (1st vs 2nd Order)
fig4 = figure('Name', 'ODE Solver - Thermal Control Accuracy Sweep', ...
              'Position', [180, 160, 950, 560], 'Color', 'w');

sig_T_mK = sigma_T_vec * 1e3;
% Highlight typical laboratory hardware accuracy regimes
patch([1, 5, 5, 1], [-65, -65, -15, -15], [0.90, 0.98, 0.90], ...
      'EdgeColor', 'none', 'DisplayName', 'Ultra-High Precision Closed Loop (PID, \sigma_T \le 5 mK)');
hold on;
patch([5, 25, 25, 5], [-65, -65, -15, -15], [0.98, 0.98, 0.85], ...
      'EdgeColor', 'none', 'DisplayName', 'Standard Commercial Lab TEC (\sigma_T \sim 5-25 mK)');
patch([25, 100, 100, 25], [-65, -65, -15, -15], [1.00, 0.90, 0.90], ...
      'EdgeColor', 'none', 'DisplayName', 'Open-Loop / Uncooled Driver (\sigma_T \ge 25 mK)');

plot(sig_T_mK, NMSE_sweep_1st_dB, 'b-o', 'LineWidth', 2.2, 'MarkerSize', 7, ...
     'MarkerFaceColor', 'b', 'DisplayName', '1st-Order ODE Solver (Single MRR)');
plot(sig_T_mK, NMSE_sweep_2nd_dB, 'r-s', 'LineWidth', 2.4, 'MarkerSize', 7, ...
     'MarkerFaceColor', 'r', 'DisplayName', '2nd-Order ODE Solver (Cascaded MRRs, 4 Heaters)');

yline(10*log10(NMSE_1st_nom), 'b:', 'LineWidth', 1.4, ...
      'DisplayName', sprintf('1st-Order Theoretical Floor (%.1f dB)', 10*log10(NMSE_1st_nom)));
yline(10*log10(NMSE_casc_nom), 'r:', 'LineWidth', 1.4, ...
      'DisplayName', sprintf('2nd-Order Theoretical Floor (%.1f dB)', 10*log10(NMSE_casc_nom)));

grid on; box on;
xlim([1, 100]); ylim([-65, -15]);
set(gca, 'XScale', 'log');
xlabel('Heater Temperature Control Accuracy \sigma_T [mK] (Log Scale)');
ylabel('Average Through-Port NMSE [dB]');
title('Thermal Precision Sensitivity: 1st-Order vs. 2nd-Order Cascaded ODE Solver');
legend('Location', 'southeast', 'FontSize', 8.5);

fprintf('========================================================================================\n');
fprintf('Simulation complete. Figures 1, 2, 3, and 4 successfully generated.\n');
fprintf('========================================================================================\n');
