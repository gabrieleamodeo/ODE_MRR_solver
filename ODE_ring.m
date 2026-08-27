%% Optical Ordinary Differential Equation (ODE) Solver via Silicon Photonic MRR
% Based on: "Compact tunable silicon photonic ODE solver", Optics Express (2014)
% Implementation of the Through-Port First-Order All-Optical ODE Solver

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

% Initial power coupling coefficients of the two couplers
kappa1     = 0.08;               % Coupler 1 (bus/through port)
kappa2     = 0.04;               % Coupler 2 (drop port)

% Ring cavity round-trip time
tau_rt     = ng * L / c;         % [s]

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

fprintf('====================================================\n');
fprintf('        Silicon Photonic ODE Solver Parameters      \n');
fprintf('====================================================\n');
fprintf('Resonant Wavelength (lambda0) : %.3f nm\n', lambda0 * 1e9);
fprintf('Ring Circumference (L)        : %.2f um\n', L * 1e6);
fprintf('Group Index (ng)              : %.4f\n', ng);
fprintf('Round-trip Time (tau_rt)      : %.3f ps\n', tau_rt * 1e12);
fprintf('Round-trip Loss (eta)         : %.4f (%.2f dB/cm)\n', eta, alpha_dB_cm);
fprintf('Coupling Coefficients         : kappa1 = %.4f, kappa2 = %.4f\n', kappa1, kappa2);
fprintf('----------------------------------------------------\n');
fprintf('Quality Factors               : Qi = %.2e, Qe1 = %.2e, Qe2 = %.2e\n', Qi, Qe1, Qe2);
fprintf('Total Quality Factor (Q)      : %.2e\n', Q_total);
fprintf('Decay Rates (gamma)           : gamma_i  = %.3e rad/s (%.2f GHz)\n', gamma_i, gamma_i/(2*pi*1e9));
fprintf('                                gamma_e1 = %.3e rad/s (%.2f GHz)\n', gamma_e1, gamma_e1/(2*pi*1e9));
fprintf('                                gamma_e2 = %.3e rad/s (%.2f GHz)\n', gamma_e2, gamma_e2/(2*pi*1e9));
fprintf('ODE Coefficient a0            : %.4e rad/s (Bandwidth: %.2f GHz)\n', a0, a0/(pi*1e9));
fprintf('ODE Coefficient b0            : %.4e rad/s\n', b0);
if b0 < 0
    fprintf('Coupling Regime               : Over-coupled (b0 < 0)\n');
elseif b0 == 0
    fprintf('Coupling Regime               : Critically coupled (b0 = 0)\n');
else
    fprintf('Coupling Regime               : Under-coupled (b0 > 0)\n');
end
fprintf('====================================================\n\n');

%% 3. Signal Definition (10 Gb/s Gaussian-like Optical Pulse)
FWHM       = 45e-12;             % Full Width at Half Maximum: 45 ps (~40-50 ps)
t0         = 0.3e-9;             % Center of the pulse [s] (300 ps)

% Simulation time grid
t_start    = -0.5e-9;            % -0.5 ns
t_end      = 1.5e-9;             % 1.5 ns
N          = 65536;              % Number of grid points (power of 2 for FFT)
t          = linspace(t_start, t_end, N);
dt         = t(2) - t(1);

% Optical complex envelope x(t) and its analytical time derivative dx(t)/dt
% Intensity: I_in(t) = exp(-4*log(2)*((t - t0)/FWHM)^2) -> Field: x(t) = sqrt(I_in(t))
x_fun      = @(t_val) exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);
dxdt_fun   = @(t_val) (-4 * log(2) * (t_val - t0) / (FWHM^2)) .* exp(-2 * log(2) * ((t_val - t0) / FWHM).^2);

x          = x_fun(t);           % Input field envelope
dxdt       = dxdt_fun(t);        % Analytical input derivative

%% 4. Numerical Ground-Truth ODE Solution (ode45)
% Target ODE: dy/dt + a0*y = dx/dt + b0*x => dy/dt = dx/dt + b0*x(t) - a0*y(t)
odefun     = @(t_val, y_val) dxdt_fun(t_val) + b0 * x_fun(t_val) - a0 * y_val;
ode_options= odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

y0         = 0;                  % Initial condition at t_start
[~, y_ideal_col] = ode45(odefun, t, y0, ode_options);
y_ideal    = y_ideal_col.';      % Row vector matching t

%% 5. Photonic ODE Solver: Frequency-Domain Transfer Function (FFT/IFFT)
% Baseband Through-Port Transfer Function: T(w) = (1j*w + b0) / (1j*w + a0)
fs         = 1 / dt;             % Sampling frequency [Hz]
f          = linspace(-fs/2, fs/2 - fs/N, N); % Centered frequency axis [Hz]
w          = 2 * pi * f;         % Baseband angular frequency [rad/s]

T_w        = (1j * w + b0) ./ (1j * w + a0);  % MRR Baseband Through-port Transfer Function

% Optical filtering via FFT
X_w        = fftshift(fft(x));     % Shifted input spectrum
Y_w        = X_w .* T_w;           % Output spectrum after MRR through-port
y_opt      = ifft(ifftshift(Y_w)); % Output field in time domain

%% 6. Power Intensities & Normalized Mean Square Error (NMSE)
I_in       = abs(x).^2;            % Input intensity
I_ideal    = abs(y_ideal).^2;      % Ground-truth output intensity (ode45)
I_opt      = abs(y_opt).^2;        % Optical solver output intensity (MRR)

% NMSE between detected optical intensities:
NMSE       = sum((I_ideal - I_opt).^2) / sum(I_ideal.^2);
NMSE_dB    = 10 * log10(NMSE);

fprintf('====================================================\n');
fprintf('                Verification Results                \n');
fprintf('====================================================\n');
fprintf('Normalized Mean Square Error (NMSE) : %.4e\n', NMSE);
fprintf('NMSE in Percentage                  : %.6f %%\n', NMSE * 100);
fprintf('NMSE in dB                          : %.2f dB\n', NMSE_dB);
fprintf('====================================================\n\n');

%% 7. Visualization & Plotting
set(0, 'DefaultAxesFontSize', 11, 'DefaultLineLineWidth', 1.8);

% Figure 1: Time-Domain Waveforms and Solver Verification
figure('Name', 'Through-Port Photonic ODE Solver - Time Domain', 'Position', [100, 100, 900, 800], 'Color', 'w');

% Subplot 1: Input Waveform
subplot(3, 1, 1);
plot(t * 1e12, I_in, 'Color', [0.1, 0.1, 0.1], 'DisplayName', 'Input Intensity |x(t)|^2');
hold on;
plot(t * 1e12, x, '--', 'Color', [0.3, 0.5, 0.9], 'LineWidth', 1.4, 'DisplayName', 'Input Field x(t)');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Amplitude / Intensity');
title('(a) Input 10 Gb/s Gaussian Pulse (FWHM = 45 ps)');
legend('Location', 'northeast');

% Subplot 2: Output Intensity Comparison (ode45 vs MRR)
subplot(3, 1, 2);
plot(t * 1e12, I_ideal, 'k-', 'DisplayName', 'Ground Truth (ode45): |y_{ideal}(t)|^2');
hold on;
plot(t * 1e12, I_opt, 'r--', 'DisplayName', 'Optical MRR Solver: |y_{opt}(t)|^2');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Output Intensity [a.u.]');
title(sprintf('(b) Output Intensity Comparison: ode45 vs. MRR Through-Port (NMSE = %.2e)', NMSE));
legend('Location', 'northeast');

% Subplot 3: Absolute Error & Field Comparison
subplot(3, 1, 3);
% plot(t * 1e12, abs(I_ideal - I_opt), 'm-', 'DisplayName', '|I_{ideal}(t) - I_{opt}(t)|');
hold on;
plot(t * 1e12, real(y_ideal), 'k:', 'LineWidth', 1.4, 'DisplayName', 'Re\{y_{ideal}(t)\}');
% plot(t * 1e12, real(y_opt), 'b-.', 'LineWidth', 1.2, 'DisplayName', 'Re\{y_{opt}(t)\}');
grid on; box on;
xlim([t_start*1e12, 1.0e3]);
xlabel('Time [ps]');
ylabel('Field / Intensity Error');
title('(c) Intensity Error & Output Field Envelopes');
legend('Location', 'northeast');

% Figure 2: Frequency Response & Spectral Characteristics
figure('Name', 'Photonic ODE Solver - Frequency Domain', 'Position', [150, 150, 900, 600], 'Color', 'w');

f_GHz = f / 1e9;
freq_limit = 50; % Plot range +/- 50 GHz

subplot(2, 1, 1);
plot(f_GHz, 10 * log10(abs(T_w).^2), 'r-', 'DisplayName', '|T(\omega)|^2 (Through Port)');
hold on;
% Normalized input spectrum for comparison
input_spec_dB = 10 * log10(abs(X_w).^2 / max(abs(X_w).^2));
plot(f_GHz, input_spec_dB, 'k--', 'LineWidth', 1.4, 'DisplayName', 'Input Pulse Spectrum (Normalized)');
yline(10 * log10(abs(b0/a0)^2), 'b:', 'LineWidth', 1.2, 'DisplayName', sprintf('Resonance Notch (|b_0/a_0|^2 = %.2f dB)', 10*log10(abs(b0/a0)^2)));
grid on; box on;
xlim([-freq_limit, freq_limit]);
ylim([-35, 5]);
xlabel('Frequency Offset f - f_0 [GHz]');
ylabel('Magnitude [dB]');
title(sprintf('(a) Through-Port Transfer Function |T(\\omega)|^2 and Input Spectrum (a_0/2\\pi = %.2f GHz)', a0/(2*pi*1e9)));
legend('Location', 'southeast');

subplot(2, 1, 2);
plot(f_GHz, angle(T_w) * (180 / pi), 'b-', 'DisplayName', '\angle T(\omega) Phase');
grid on; box on;
xlim([-freq_limit, freq_limit]);
ylim([-190, 190]);
xlabel('Frequency Offset f - f_0 [GHz]');
ylabel('Phase [degrees]');
title('(b) Through-Port Phase Response \angle T(\omega)');
legend('Location', 'southeast');
