# Silicon Photonic Ordinary Differential Equation (ODE) Solver

A comprehensive MATLAB/Octave simulation framework for **all-optical ordinary differential equation (ODE) solvers** based on tunable silicon-on-insulator (SOI) microring resonators (MRRs). 

This repository models, simulates, and analyzes both **first-order** and cascaded **second-order** optical ODE solvers, investigating physical waveguide losses, thermal tuning boundaries, hardware thermal noise, and inter-cavity resonance mismatches.

---

## Table of Contents
1. [Theoretical Foundations](#theoretical-foundations)
   - [First-Order Optical ODE Solver](#1-first-order-optical-ode-solver)
   - [MZI Tunable Directional Couplers](#2-mzi-tunable-directional-couplers)
   - [Cascaded Second-Order Optical ODE Solver](#3-cascaded-second-order-optical-ode-solver)
2. [Repository Scripts & Modules](#repository-scripts--modules)
   - [`ODE_ring.m`](#1-ode_ringm---baseline-first-order-solver)
   - [`ODE_ring_cavity_loss_analysis.m`](#2-ode_ring_cavity_loss_analysism---waveguide-loss-sweep)
   - [`ODE_ring_thermal_tuning_bounds.m`](#3-ode_ring_thermal_tuning_boundsm---reachable-ode-bounds)
   - [`ODE_ring_thermal_noise.m`](#4-ode_ring_thermal_noisem---1st-order-thermal-noise--monte-carlo)
   - [`ODE_ring_second_order_thermal_noise.m`](#5-ode_ring_second_order_thermal_noisem---2nd-order-cascade--independent-noise)
3. [Requirements & Execution](#requirements--execution)
4. [References](#references)

---

## Theoretical Foundations

Based on **J. Wu et al., *"Compact tunable silicon photonic differential-equation solver for general linear time-invariant systems,"* Optics Express 22(21), 26254–26264 (2014)**.

### 1. First-Order Optical ODE Solver

An add-drop microring resonator (MRR) excited at carrier angular frequency $\omega_0$ performs analog computation on the complex envelope $x(t)$ of an input optical pulse. The through-port optical output $y(t)$ satisfies the first-order linear time-invariant (LTI) differential equation:

$$\frac{dy(t)}{dt} + a_0 y(t) = \frac{dx(t)}{dt} + b_0 x(t)$$

Taking the Fourier transform yields the baseband optical transfer function:

$$T_1(\omega) = \frac{j\omega + b_0}{j\omega + a_0}$$

From Coupled-Mode Theory (CMT), the ODE coefficients $a_0$ and $b_0$ are directly mapped to the physical decay rates of the cavity:

$$a_0 = \gamma_i + \gamma_{e1} + \gamma_{e2}, \qquad b_0 = \gamma_i + \gamma_{e2} - \gamma_{e1}$$

where:
* $\gamma_i = \frac{c \alpha_{\text{lin}}}{2 n_g}$ is the intrinsic round-trip power loss rate [rad/s], with $\alpha_{\text{lin}} = \alpha_{\text{dB/m}} \frac{\ln(10)}{10}$.
* $\gamma_{e1} = -\frac{c}{2 n_g L} \ln(1 - \kappa_1)$ is the external coupling rate to the input/through bus waveguide [rad/s].
* $\gamma_{e2} = -\frac{c}{2 n_g L} \ln(1 - \kappa_2)$ is the external coupling rate to the drop bus waveguide [rad/s].
* $L = 178.98\ \mu\text{m}$ is the cavity circumference, $n_g = 4.1850$ is the waveguide group index, and $\lambda_0 = 1550.391\text{ nm}$.

### 2. MZI Tunable Directional Couplers

The power coupling ratios $\kappa_1$ and $\kappa_2$ are dynamically tuned via asymmetric Mach-Zehnder Interferometers (MZIs) integrated with microheaters:

$$\kappa(\Delta\phi) = 2 \kappa_0 (1 - \kappa_0) (1 + \cos\Delta\phi)$$

* $\kappa_0 = 0.0441$ is the power coupling ratio of the internal directional coupler (FDTD design).
* $\kappa_{\max} = 4 \kappa_0 (1 - \kappa_0) \approx 0.168621$ is the maximum achievable coupling.
* The phase shift $\Delta\phi$ is tuned thermo-optically:
  $$\Delta\phi(T) = \Delta\phi_0 + \frac{2\pi}{\lambda_0} \frac{dn}{dT} L_b \Delta T$$
  with $\frac{dn}{dT} = 1.86 \times 10^{-4}\text{ K}^{-1}$ (thermo-optic coefficient of silicon) and $L_b$ the heater length.


## Repository Scripts & Modules

| Script | Purpose & Description | Key Output Figures |
| :--- | :--- | :--- |
| [`ODE_ring.m`](ODE_ring.m) | **Baseline Verification:** Compares the Ideal Coupled-Mode Theory (CMT) ODE solver against the exact physical add-drop MRR model with a 10 Gb/s Gaussian pulse (FWHM = 45 ps). | Time-domain waveforms, Baseband spectrum, Periodic FSR spectrum. |
| [`ODE_ring_cavity_loss_analysis.m`](ODE_ring_cavity_loss_analysis.m) | **Waveguide Loss Analysis:** Sweeps propagation loss $\alpha \in [2, 20]\text{ dB/cm}$. Explains the critical coupling condition at $\alpha \approx 10.33\text{ dB/cm}$ and NMSE behavior. | Cavity rates vs. $\alpha$, Realized $a_0, b_0$, Through-port NMSE curve, Time-domain waveforms at distinct loss regimes. |
| [`ODE_ring_thermal_tuning_bounds.m`](ODE_ring_thermal_tuning_bounds.m) | **Reachable ODE Coefficient Bounds:** Determines the physical $[a_0, b_0]$ tuning bounds as a function of microheater temperature sweeps $\Delta T \in [0, 100]\text{ K}$. | 2D parameter boundary map $[a_0, b_0]$ in $\text{rad/s}$, Coupler cross-coupling $\kappa(T)$, Temperature scale maps. |
| [`ODE_ring_thermal_noise.m`](ODE_ring_thermal_noise.m) | **First-Order Thermal Noise & Monte Carlo:** Simulates DAC noise, TEC thermal drift, and cross-talk ($\sigma_T \in [1, 100]\text{ mK}$) on a single MRR. | Realization waveforms, Baseband detuned spectrum, $N=200$ Monte Carlo confidence bands, $\sigma_T$ sweep. |
| [`ODE_ring_second_order_thermal_noise.m`](ODE_ring_second_order_thermal_noise.m) | **Second-Order Cascaded ODE Solver:** Simulates two cascaded MRRs with 4 independent heaters and 2 independent cavity drifts. Compares 1st vs. 2nd-order performance. | Waveform confidence bands, Asymmetric notch splitting, NMSE statistical histograms, Comparative thermal sweep. |

---

## Requirements & Execution

### Prerequisites
* **MATLAB** (R2018b or later) or **GNU Octave** (v6.0+ with `signal` package).
* No proprietary toolboxes are required; all FFTs, Coupled-Mode Theory models, and transfer function matrix calculations are self-contained.

### Running Simulations
Clone the repository:
```bash
git clone https://github.com/gabrieleamodeo/ODE_MRR_solver.git
cd ODE_MRR_solver
```

In MATLAB or Octave:
```matlab
% 1. Run baseline verification (First-Order ODE vs. Physical MRR)
run('ODE_ring.m')

% 2. Run waveguide loss analysis (Sweep alpha in [2, 20] dB/cm)
run('ODE_ring_cavity_loss_analysis.m')

% 3. Run thermal tuning reachable coefficient map
run('ODE_ring_thermal_tuning_bounds.m')

% 4. Run single-ring thermal noise Monte Carlo simulation
run('ODE_ring_thermal_noise.m')

% 5. Run second-order cascaded ODE solver with 4 independent heaters
run('ODE_ring_second_order_thermal_noise.m')
```

---

## Documentation & Reports

* **Reference Paper:** [`docs/4.6. Compact tunable silicon photonic ODE solver.pdf`](docs/4.6.%20Compact%20tunable%20silicon%20photonic%20ODE%20solver.pdf) — J. Wu et al. (Optics Express 2014).
* **Technical Report:** [`docs/Report_ODE_MRR_Solver.pdf`](docs/Report_ODE_MRR_Solver.pdf) — Comprehensive technical report detailing the design, derivations, thermal noise models, and cascaded second-order ODE analysis.

---

## References

1. **J. Wu, P. Cao, X. Hu, X. Jiang, T. Pan, Y. Yang, C. Qiu, C. Tremblay, and Y. Su**, *"Compact tunable silicon photonic differential-equation solver for general linear time-invariant systems,"* **Optics Express**, vol. 22, no. 21, pp. 26254–26264, 2014. [DOI: 10.1364/OE.22.026254](https://doi.org/10.1364/OE.22.026254).
2. **F. Liu, T. Wang, L. Qiang, T. Ye, Z. Zhang, M. Qiu, and Y. Su**, *"Compact optical temporal differentiator based on silicon microring resonator,"* **Optics Express**, vol. 16, no. 20, pp. 15880–15886, 2008.
3. **B. E. A. Saleh and M. C. Teich**, *Fundamentals of Photonics*, 3rd ed. Hoboken, NJ: Wiley, 2019.
