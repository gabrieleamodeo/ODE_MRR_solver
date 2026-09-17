%% Silicon Photonic MRR ODE Solver: Thermal Tuning Range & Coefficient Extrema
% Reference: J. Wu et al., "Compact tunable silicon photonic differential-equation
%            solver for general linear time-invariant systems," Optics Express 22(21),
%            26254-26264 (2014).
%
% Theoretical & Physical Purpose:
% -------------------------------
% Given the typical silicon waveguide propagation loss of the paper (alpha = 8.0 dB/cm),
% this script computes the exact minimum and maximum achievable values of the
% first-order ODE coefficients [a0, b0] as a function of the minimum and maximum
% temperatures (T_min, T_max) applied to the thermo-optic microheaters controlling
% the two interferometric couplers.
%
% Governing Equations (Coupled-Mode Theory):
%   dy(t)/dt + a0*y(t) = dx(t)/dt + b0*x(t)
%
%   a0(T1, T2) = gamma_i + gamma_e1(T1) + gamma_e2(T2)
%   b0(T1, T2) = gamma_i + gamma_e2(T2) - gamma_e1(T1)
%
% where:
%   gamma_i       = (c * alpha_lin) / (2 * ng)  [Constant, fixed by alpha = 8 dB/cm]
%   gamma_e1(T1)  = - (c * ln(1 - kappa1(T1))) / (2 * ng * L)
%   gamma_e2(T2)  = - (c * ln(1 - kappa2(T2))) / (2 * ng * L)
%   kappa1,2(T)   = 2 * kappa0 * (1 - kappa0) * [1 + cos(Delta_phi1,2(T))]
%   Delta_phi1(T) = Delta_phi1_0 + (2*pi/lambda0) * (dn/dT) * Lb1 * (T1 - T_ambient)
%   Delta_phi2(T) = Delta_phi2_0 + (2*pi/lambda0) * (dn/dT) * Lb2 * (T2 - T_ambient)

clear; close all; clc;

%% 1. Baseline Physical Constants & Waveguide Parameters (Wu et al. 2014)
c             = 2.99792458e8;       % Speed of light in vacuum [m/s]
lambda0       = 1550.391e-9;        % Operating optical wavelength [m] (1550.391 nm)
f0            = c / lambda0;        % Optical carrier frequency [Hz] (~193.366 THz)
w0            = 2 * pi * f0;        % Optical carrier angular frequency [rad/s]
L             = 178.98e-6;          % Ring circumference [m] (178.98 um)
ng            = 4.1850;             % Silicon waveguide group index (TE mode)
tau_rt        = ng * L / c;         % Round-trip time [s] (~2.498 ps)
FSR           = 1 / tau_rt;         % Free Spectral Range [Hz] (~400.24 GHz)

% Waveguide Propagation Loss from Paper
alpha_dB_cm   = 8.0;                % Waveguide power loss factor [dB/cm]
alpha_dB_m    = alpha_dB_cm * 100;  % [dB/m]
alpha_lin     = alpha_dB_m * log(10) / 10; % Linear attenuation coefficient [1/m] (~184.2 1/m)
gamma_i       = c * alpha_lin / (2 * ng);  % Intrinsic cavity decay rate [rad/s] (~6.598e9 rad/s)
eta_roundtrip = 1 - exp(-alpha_lin * L);    % Round-trip intrinsic power loss (~3.25%)
Qi            = w0 / (2 * gamma_i);        % Intrinsic cavity quality factor (~9.21e4)

% MZI Tunable Interferometric Coupler Parameters (Section 3 of Paper)
kappa0        = 0.0441;             % Directional coupler splitting power ratio (FDTD)
kappa_max     = 4 * kappa0 * (1 - kappa0); % Maximum MZI power coupling (~0.168621)
gamma_e_max   = - c * log(1 - kappa_max) / (2 * ng * L); % Max coupling rate (~3.696e10 rad/s)

% Thermo-Optic Properties & Heater Arm Lengths (Wu et al. 2014)
dn_dT         = 1.86e-4;            % Thermo-optic coefficient of silicon [1/K]
Lb1           = 116.80e-6;          % Coupler 1 bus arm heater length [m] (116.80 um)
Lr1           = 47.12e-6;           % Coupler 1 ring arm length [m] (47.12 um)
Lb2           = 47.12e-6;           % Coupler 2 bus arm heater length [m] (47.12 um)
Lr2           = 47.12e-6;           % Coupler 2 ring arm length [m] (47.12 um)

% Geometrical Path Length Differences
dL1           = Lb1 - Lr1;          % Arm length difference coupler 1: 69.68 um
dL2           = Lb2 - Lr2;          % Arm length difference coupler 2: 0.00 um (symmetric arms)

% Static Geometrical Phase Offsets at Ambient Temperature (dT = 0)
dphi1_0       = (2 * pi / lambda0) * ng * dL1; % Static phase offset coupler 1 [rad]
dphi2_0       = (2 * pi / lambda0) * ng * dL2; % Static phase offset coupler 2 = 0 [rad]

% Characteristic Temperature Tuning Scales
% Temperature increase required for pi phase shift (minimum to maximum coupling):
dT_pi_1       = lambda0 / (2 * dn_dT * Lb1); % ~35.68 K for Heater 1
dT_pi_2       = lambda0 / (2 * dn_dT * Lb2); % ~88.45 K for Heater 2
% Temperature increase required for full 2*pi phase cycle:
dT_2pi_1      = lambda0 / (dn_dT * Lb1);     % ~71.37 K for Heater 1
dT_2pi_2      = lambda0 / (dn_dT * Lb2);     % ~176.90 K for Heater 2

%% 2. User-Defined Operating Temperature Window
% You can adjust these temperatures to evaluate your specific lab setup!
T_ambient_C   = 20.0;               % Ambient / baseline chip temperature [deg C]
T_min_C       = 20.0;               % Minimum achievable heater temperature [deg C]
T_max_C       = 120.0;              % Maximum allowable heater temperature [deg C]

% Temperature elevations relative to ambient [K]
dT_min        = max(0, T_min_C - T_ambient_C);
dT_max        = max(dT_min, T_max_C - T_ambient_C);

%% 3. Analytical Global Physical Extrema (Device Theoretical Boundaries)
% Regardless of temperature, the physical interferometric couplers can only
% span power couplings in the range kappa in [0, kappa_max].
% Because gamma_e(kappa) is strictly monotonic with kappa:
%   gamma_e in [0, gamma_e_max]
%
% Absolute Global Extrema for alpha = 8.0 dB/cm:
% a0 = gamma_i + gamma_e1 + gamma_e2:
a0_global_min = gamma_i;                                % kappa1 = 0, kappa2 = 0
a0_global_max = gamma_i + 2 * gamma_e_max;              % kappa1 = kappa_max, kappa2 = kappa_max

% b0 = gamma_i + gamma_e2 - gamma_e1:
b0_global_min = gamma_i - gamma_e_max;                  % kappa1 = kappa_max, kappa2 = 0 (over-coupled)
b0_global_max = gamma_i + gamma_e_max;                  % kappa1 = 0, kappa2 = kappa_max (under-coupled)

% The 4 Vertices of the Complete Achievable Parallelogram in (a0, b0) space:
V1 = [a0_global_min,                      gamma_i];                   % (k1=0, k2=0)
V2 = [gamma_i + gamma_e_max,              gamma_i - gamma_e_max];     % (k1=max, k2=0)
V3 = [a0_global_max,                      gamma_i];                   % (k1=max, k2=max)
V4 = [gamma_i + gamma_e_max,              gamma_i + gamma_e_max];     % (k1=0, k2=max)

%% 4. Achievable Extrema for the Specific Temperature Window [T_min, T_max]
% Sample the [dT_min, dT_max] temperature window finely
N_pts_T       = 300;
dT1_vec       = linspace(dT_min, dT_max, N_pts_T);
dT2_vec       = linspace(dT_min, dT_max, N_pts_T);

% Coupler phase shifts across the operating temperature span
phi1_vec      = dphi1_0 + (2 * pi / lambda0) * dn_dT * Lb1 * dT1_vec;
phi2_vec      = dphi2_0 + (2 * pi / lambda0) * dn_dT * Lb2 * dT2_vec;

% Coupler power coupling ranges achievable within [T_min, T_max]
k1_T_vec      = 2 * kappa0 * (1 - kappa0) * (1 + cos(phi1_vec));
k2_T_vec      = 2 * kappa0 * (1 - kappa0) * (1 + cos(phi2_vec));

k1_win_min    = min(k1_T_vec);
k1_win_max    = max(k1_T_vec);
k2_win_min    = min(k2_T_vec);
k2_win_max    = max(k2_T_vec);

% Corresponding decay rate ranges
ge1_win_min   = - c * log(max(1 - k1_win_min, 1e-12)) / (2 * ng * L);
ge1_win_max   = - c * log(max(1 - k1_win_max, 1e-12)) / (2 * ng * L);
ge2_win_min   = - c * log(max(1 - k2_win_min, 1e-12)) / (2 * ng * L);
ge2_win_max   = - c * log(max(1 - k2_win_max, 1e-12)) / (2 * ng * L);

% Achievable Extrema within [T_min, T_max]:
a0_win_min    = gamma_i + ge1_win_min + ge2_win_min;
a0_win_max    = gamma_i + ge1_win_max + ge2_win_max;

b0_win_min    = gamma_i + ge2_win_min - ge1_win_max;
b0_win_max    = gamma_i + ge2_win_max - ge1_win_min;

%% 5. Parametric Sweep: Extrema as a Continuous Function of Delta_T_max
% Here we sweep the upper thermal tuning limit Delta_T_max from 0 to 150 K
% to show how the bounds expand until saturation.
dT_sweep_max  = linspace(0, 150, 151); % 0 to 150 K temperature elevation
N_sweep       = length(dT_sweep_max);

a0_min_sweep  = zeros(1, N_sweep);
a0_max_sweep  = zeros(1, N_sweep);
b0_min_sweep  = zeros(1, N_sweep);
b0_max_sweep  = zeros(1, N_sweep);
k1_min_sweep  = zeros(1, N_sweep);
k1_max_sweep  = zeros(1, N_sweep);
k2_min_sweep  = zeros(1, N_sweep);
k2_max_sweep  = zeros(1, N_sweep);

for s = 1:N_sweep
    dt_lim = dT_sweep_max(s);
    if dt_lim == 0
        p1 = dphi1_0;
        p2 = dphi2_0;
    else
        dt_sub = linspace(0, dt_lim, 100);
        p1 = dphi1_0 + (2 * pi / lambda0) * dn_dT * Lb1 * dt_sub;
        p2 = dphi2_0 + (2 * pi / lambda0) * dn_dT * Lb2 * dt_sub;
    end
    k1_sub = 2 * kappa0 * (1 - kappa0) * (1 + cos(p1));
    k2_sub = 2 * kappa0 * (1 - kappa0) * (1 + cos(p2));
    
    k1_min_sweep(s) = min(k1_sub);
    k1_max_sweep(s) = max(k1_sub);
    k2_min_sweep(s) = min(k2_sub);
    k2_max_sweep(s) = max(k2_sub);
    
    ge1_min_s = - c * log(max(1 - k1_min_sweep(s), 1e-12)) / (2 * ng * L);
    ge1_max_s = - c * log(max(1 - k1_max_sweep(s), 1e-12)) / (2 * ng * L);
    ge2_min_s = - c * log(max(1 - k2_min_sweep(s), 1e-12)) / (2 * ng * L);
    ge2_max_s = - c * log(max(1 - k2_max_sweep(s), 1e-12)) / (2 * ng * L);
    
    a0_min_sweep(s) = gamma_i + ge1_min_s + ge2_min_s;
    a0_max_sweep(s) = gamma_i + ge1_max_s + ge2_max_s;
    b0_min_sweep(s) = gamma_i + ge2_min_s - ge1_max_s;
    b0_max_sweep(s) = gamma_i + ge2_max_s - ge1_min_s;
end

%% 6. 2D Meshgrid Parameter Space across [T_min, T_max]
[DT1, DT2]    = meshgrid(dT1_vec, dT2_vec);
Phi1_2D       = dphi1_0 + (2 * pi / lambda0) * dn_dT * Lb1 * DT1;
Phi2_2D       = dphi2_0 + (2 * pi / lambda0) * dn_dT * Lb2 * DT2;

K1_2D         = 2 * kappa0 * (1 - kappa0) * (1 + cos(Phi1_2D));
K2_2D         = 2 * kappa0 * (1 - kappa0) * (1 + cos(Phi2_2D));

Ge1_2D        = - c * log(max(1 - K1_2D, 1e-12)) / (2 * ng * L);
Ge2_2D        = - c * log(max(1 - K2_2D, 1e-12)) / (2 * ng * L);

A0_2D         = gamma_i + Ge1_2D + Ge2_2D;
B0_2D         = gamma_i + Ge2_2D - Ge1_2D;

%% 7. Formatted Console Summary & Physical Insights
fprintf('========================================================================================\n');
fprintf('       SILICON PHOTONIC MRR ODE SOLVER: THERMAL TUNING & COEFFICIENT BOUNDS             \n');
fprintf('       Reference: Wu et al., Optics Express 22(21), 26254-26264 (2014)                  \n');
fprintf('========================================================================================\n');
fprintf('Waveguide Propagation Loss (alpha)  : %.2f dB/cm (Linear: %.2f m^-1)\n', alpha_dB_cm, alpha_lin);
fprintf('Cavity Round-Trip Loss (eta)        : %.3f %% (Intrinsic Qi = %.2e)\n', eta_roundtrip * 100, Qi);
fprintf('Fixed Intrinsic Decay Rate (gamma_i): %.4e rad/s\n', gamma_i);
fprintf('Max MZI Power Coupling (kappa_max)  : %.6f (Directional coupler kappa0 = %.4f)\n', kappa_max, kappa0);
fprintf('Max Coupling Decay Rate (gamma_max) : %.4e rad/s\n', gamma_e_max);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('THERMO-OPTIC HEATER CHARACTERISTICS:\n');
fprintf('  Heater 1 (Through Coupler Bus Arm): Lb1 = %.2f um, Lr1 = %.2f um (Delta_L1 = %.2f um)\n', ...
        Lb1 * 1e6, Lr1 * 1e6, dL1 * 1e6);
fprintf('    -> Delta_T for pi phase shift   : %.2f K (%.1f deg C above ambient)\n', dT_pi_1, T_ambient_C + dT_pi_1);
fprintf('    -> Delta_T for 2*pi phase shift : %.2f K (%.1f deg C above ambient)\n', dT_2pi_1, T_ambient_C + dT_2pi_1);
fprintf('  Heater 2 (Drop Coupler Bus Arm)   : Lb2 = %.2f um, Lr2 = %.2f um (Delta_L2 = 0.00 um)\n', ...
        Lb2 * 1e6, Lr2 * 1e6);
fprintf('    -> Delta_T for pi phase shift   : %.2f K (%.1f deg C above ambient)\n', dT_pi_2, T_ambient_C + dT_pi_2);
fprintf('    -> Delta_T for 2*pi phase shift : %.2f K (%.1f deg C above ambient)\n', dT_2pi_2, T_ambient_C + dT_2pi_2);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('ABSOLUTE GLOBAL PHYSICAL EXTREMA (Device Architecture Limits at alpha = 8.0 dB/cm):\n');
fprintf('  Coefficient a0 (Bandwidth / Damping) :\n');
fprintf('    - Minimum a0 [kappa1=0, kappa2=0] : %.4e rad/s\n', a0_global_min);
fprintf('    - Maximum a0 [kappa1=max, kappa2=max]: %.4e rad/s\n', a0_global_max);
fprintf('    - Full a0 Tuning Span (Delta a0)  : %.4e rad/s (Dynamic Ratio = %.1fx)\n', ...
        a0_global_max - a0_global_min, a0_global_max/a0_global_min);
fprintf('  Coefficient b0 (Zero Location / Coupling Symmetry) :\n');
fprintf('    - Minimum b0 [kappa1=max, kappa2=0] : %+.4e rad/s -> Deep Over-coupled\n', ...
        b0_global_min);
fprintf('    - Maximum b0 [kappa1=0, kappa2=max] : %+.4e rad/s -> Deep Under-coupled\n', ...
        b0_global_max);
fprintf('    - Critical Coupling Zero-Crossing  : b0 = 0 rad/s (gamma_e1 - gamma_e2 = gamma_i = %.4e rad/s)\n', ...
        gamma_i);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('CONFIGURED OPERATING WINDOW: T in [%.1f, %.1f] deg C (Delta_T in [%.1f, %.1f] K):\n', ...
        T_min_C, T_max_C, dT_min, dT_max);
fprintf('  Coupler 1 Achievable Coupling     : kappa1 in [%.4f, %.4f]\n', k1_win_min, k1_win_max);
fprintf('  Coupler 2 Achievable Coupling     : kappa2 in [%.4f, %.4f]\n', k2_win_min, k2_win_max);
fprintf('  Achieved a0 Range                 : [%.4e, %.4e] rad/s\n', a0_win_min, a0_win_max);
fprintf('  Achieved b0 Range                 : [%+.4e, %+.4e] rad/s\n', b0_win_min, b0_win_max);
if b0_win_min < 0 && b0_win_max > 0
    fprintf('  Coupling Regimes Accessible       : BOTH Over-coupled (b0 < 0) and Under-coupled (b0 > 0)!\n');
    fprintf('                                      Includes Critical Coupling Differentiator Mode (b0 = 0).\n');
elseif b0_win_max <= 0
    fprintf('  Coupling Regimes Accessible       : Exclusively Over-coupled (b0 <= 0).\n');
else
    fprintf('  Coupling Regimes Accessible       : Exclusively Under-coupled (b0 >= 0).\n');
end
fprintf('========================================================================================\n\n');

%% 8. Visualization: Figure 1 - Evolution of Extrema vs. Maximum Temperature
set(0, 'DefaultAxesFontSize', 10.5, 'DefaultLineLineWidth', 1.8);

fig1 = figure('Name', 'MRR ODE Solver - Thermal Tuning Range & Extrema', ...
              'Position', [60, 40, 950, 850], 'Color', 'w');

% (a) Power Coupling Coefficients vs Temperature Elevation
subplot(3, 1, 1);
plot(dT_sweep_max, k1_T_vec(1:N_sweep), 'b-', 'LineWidth', 2.2, 'DisplayName', '\kappa_1(\Delta T_1) (Heater 1, L_{b1}=116.8 \mum)');
hold on;
plot(dT_sweep_max, k2_T_vec(1:N_sweep), 'r-', 'LineWidth', 2.2, 'DisplayName', '\kappa_2(\Delta T_2) (Heater 2, L_{b2}=47.1 \mum)');
yline(kappa_max, 'k:', 'LineWidth', 1.5, 'DisplayName', sprintf('\\kappa_{max} = %.4f', kappa_max));
yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', '\kappa_{min} = 0');
xline(dT_pi_1, 'b--', 'LineWidth', 1.2, 'DisplayName', sprintf('\\Delta T_{\\pi,1} = %.1f K', dT_pi_1));
xline(dT_pi_2, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('\\Delta T_{\\pi,2} = %.1f K', dT_pi_2));
grid on; box on;
xlim([0, 150]); ylim([-0.01, 0.185]);
xlabel('Temperature Elevation \Delta T [K]');
ylabel('Power Coupling Ratio \kappa');
title('(a) Coupler Power Couplings vs. Heater Temperature');
legend('Location', 'northeast', 'FontSize', 8.5);

% (b) Min and Max a0 vs Temperature Tuning Range Delta_T_max
subplot(3, 1, 2);
a0_min_rad = a0_min_sweep;
a0_max_rad = a0_max_sweep;
fill([dT_sweep_max, fliplr(dT_sweep_max)], ...
     [a0_max_rad, fliplr(a0_min_rad)], ...
     [0.85, 0.92, 1.0], 'EdgeColor', 'none', 'DisplayName', 'Achievable a_0 Bandwidth Range');
hold on;
plot(dT_sweep_max, a0_max_rad, 'b-', 'LineWidth', 2.2, 'DisplayName', 'a_{0,max}(\Delta T_{max})');
plot(dT_sweep_max, a0_min_rad, 'b--', 'LineWidth', 2.0, 'DisplayName', 'a_{0,min}(\Delta T_{max})');
yline(a0_global_max, 'k:', 'LineWidth', 1.5, 'DisplayName', sprintf('Global a_{0,max} = %.2e rad/s', a0_global_max));
yline(a0_global_min, 'k-.', 'LineWidth', 1.5, 'DisplayName', sprintf('Global a_{0,min} = \\gamma_i = %.2e rad/s', a0_global_min));
xline(dT_max, 'm-', 'LineWidth', 1.6, 'DisplayName', sprintf('User Window \\Delta T_{max} = %.0f K', dT_max));
grid on; box on;
xlim([0, 150]); ylim([0, 9e10]);
xlabel('Maximum Temperature Tuning Range \Delta T_{max} [K]');
ylabel('ODE Coefficient a_0 [rad/s]');
title('(b) Bandwidth / Damping Coefficient a_0 Tuning Range');
legend('Location', 'east', 'FontSize', 8.5);

% (c) Min and Max b0 vs Temperature Tuning Range Delta_T_max
subplot(3, 1, 3);
b0_min_rad = b0_min_sweep;
b0_max_rad = b0_max_sweep;
fill([dT_sweep_max, fliplr(dT_sweep_max)], ...
     [b0_max_rad, fliplr(b0_min_rad)], ...
     [1.0, 0.88, 0.88], 'EdgeColor', 'none', 'DisplayName', 'Achievable b_0 Tuning Range');
hold on;
plot(dT_sweep_max, b0_max_rad, 'r-', 'LineWidth', 2.2, 'DisplayName', 'b_{0,max}(\Delta T_{max}) (Under-coupled)');
plot(dT_sweep_max, b0_min_rad, 'r--', 'LineWidth', 2.0, 'DisplayName', 'b_{0,min}(\Delta T_{max}) (Over-coupled)');
yline(0, 'k-', 'LineWidth', 1.4, 'DisplayName', 'Critical Coupling Line (b_0 = 0)');
yline(b0_global_max, 'k:', 'LineWidth', 1.5, 'DisplayName', sprintf('Global b_{0,max} = %+.2e rad/s', b0_global_max));
yline(b0_global_min, 'k-.', 'LineWidth', 1.5, 'DisplayName', sprintf('Global b_{0,min} = %.2e rad/s', b0_global_min));
xline(dT_max, 'm-', 'LineWidth', 1.6, 'DisplayName', sprintf('User Window \\Delta T_{max} = %.0f K', dT_max));
grid on; box on;
xlim([0, 150]); ylim([-3.5e10, 5e10]);
xlabel('Maximum Temperature Tuning Range \Delta T_{max} [K]');
ylabel('ODE Coefficient b_0 [rad/s]');
title('(c) Zero Location / Coupling Symmetry b_0 Tuning Range');
legend('Location', 'east', 'FontSize', 8.5);

%% 9. Visualization: Figure 2 - 2D Thermal Landscapes & Reachable (a0, b0) Space
fig2 = figure('Name', 'MRR ODE Solver - 2D Operational Space & Landscapes', ...
              'Position', [120, 100, 1100, 820], 'Color', 'w');

T1_grid_C = T_ambient_C + DT1;
T2_grid_C = T_ambient_C + DT2;

% (a) 2D Contour of a0 across the (T1, T2) Heater Control Plane
subplot(2, 2, 1);
contourf(T1_grid_C, T2_grid_C, A0_2D, 25, 'LineColor', 'none');
colorbar; colormap(gca, 'jet');
hold on;
xlabel('Heater 1 Temperature T_1 [^\circC]');
ylabel('Heater 2 Temperature T_2 [^\circC]');
title('(a) Coefficient a_0 [rad/s] vs. (T_1, T_2)');
grid on; box on;

% (b) 2D Contour of b0 across the (T1, T2) Heater Control Plane
subplot(2, 2, 2);
contourf(T1_grid_C, T2_grid_C, B0_2D, 25, 'LineColor', 'none');
colorbar; colormap(gca, 'parula');
hold on;
% Critical coupling differentiator contour (b0 = 0)
[C_crit, h_crit] = contour(T1_grid_C, T2_grid_C, B0_2D, [0 0], 'w-', 'LineWidth', 2.4);
xlabel('Heater 1 Temperature T_1 [^\circC]');
ylabel('Heater 2 Temperature T_2 [^\circC]');
title('(b) Coefficient b_0 [rad/s] & Differentiator Line (b_0 = 0)');
grid on; box on;

% (c) 2D Reachable ODE Space (a0, b0) - Full Physical Parallelogram vs Active Window
subplot(2, 2, [3, 4]);

% 1. Plot Global Theoretical Parallelogram
poly_a0 = [V1(1), V2(1), V3(1), V4(1), V1(1)];
poly_b0 = [V1(2), V2(2), V3(2), V4(2), V1(2)];
fill(poly_a0, poly_b0, [0.93, 0.93, 0.93], 'EdgeColor', [0.4 0.4 0.4], ...
     'LineWidth', 1.6, 'DisplayName', 'Global Physical Architecture Domain (\alpha = 8 dB/cm)');
hold on;

% 2. Scatter plot of points reachable within current [T_min, T_max]
scatter(A0_2D(:), B0_2D(:), 9, [0.2 0.6 0.9], 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', sprintf('Reachable in [%.0f, %.0f] ^\\circC Window', T_min_C, T_max_C));

% 3. Critical Coupling Differentiator Axis (b0 = 0)
yline(0, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Critical Coupling Axis (b_0 = 0, Differentiator)');

% 4. Shaded Over-coupled and Under-coupled zones text
text(1.2e10, -2.0e10, 'OVER-COUPLED REGIME (b_0 < 0)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0.1 0.1]);
text(1.2e10, +3.0e10, 'UNDER-COUPLED REGIME (b_0 > 0)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.1 0.3 0.8]);

% 5. Highlight Experimental System Testing Points from Wu et al. (2014)
% Testing I (unheated P1=0, P2=0): a0 = 8.0127e10 rad/s, b0 = 8.1233e9 rad/s
plot(8.0127e10, 8.1233e9, 'p', 'MarkerSize', 13, ...
     'MarkerFaceColor', [1.0 0.8 0.0], 'MarkerEdgeColor', 'k', ...
     'DisplayName', 'Paper Exp Testing I (Unheated: a_0=8.01\times10^{10}, b_0=+8.12\times10^9 rad/s)');

% Testing Set 3 (constant a0 = 6.3452e10 rad/s, varied b0)
a0_exp3 = 6.3452e10;
b0_exp3 = [-8.5515e9, -1.6055e9, 5.4834e9];
plot(a0_exp3 * ones(size(b0_exp3)), b0_exp3, 'rs', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'Paper Exp Set 3 (Constant a_0 = 6.35\times10^{10} rad/s, Varied b_0)');

% Testing Set 4 (constant b0 = 6.7874e9 rad/s, varied a0)
b0_exp4 = 6.7874e9;
a0_exp4 = [2.8127e10, 5.9454e10, 7.8791e10];
plot(a0_exp4, b0_exp4 * ones(size(a0_exp4)), 'g^', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'g', 'DisplayName', 'Paper Exp Set 4 (Constant b_0 = 6.79\times10^9 rad/s, Varied a_0)');

grid on; box on;
xlim([0, 9.0e10]); ylim([-3.5e10, 5.0e10]);
xlabel('ODE Bandwidth / Damping Coefficient a_0 [rad/s]');
ylabel('ODE Zero Location Coefficient b_0 [rad/s]');
title('(c) 2D Achievable ODE Coefficient Space (a_0, b_0) & Experimental Verification');
legend('Location', 'northwest', 'FontSize', 9.0);

fprintf('========================================================================================\n');
fprintf('Simulation complete. Figures 1 and 2 successfully generated.\n');
fprintf('========================================================================================\n');
