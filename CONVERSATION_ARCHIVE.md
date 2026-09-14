# Photonic ODE Solver: Thermal Noise Modeling & Monte Carlo Analysis
**Conversation Archive & Technical Reference**  
*Date: September 14, 2026*  
*Repository: `ODE_MRR_solver`*  
*Reference: J. Wu et al., "Compact tunable silicon photonic differential-equation solver for general linear time-invariant systems," Optics Express 22(21), 26254–26264 (2014).*

---

## 1. Overview & Session Objectives

This session addressed the design, physical modeling, implementation, and statistical analysis of **thermal control errors** in an on-chip silicon photonic microring resonator (MRR) ordinary differential equation (ODE) solver.

Key topics covered:
1. **AI Model Selection Configuration** in the Antigravity UI.
2. **Physical Modeling of Thermal Inaccuracies:** Formulating realistic temperature fluctuations in on-chip thermo-optic phase shifters (interferometric couplers and ring cavity detuning) with typical laboratory hardware control limits ($\sigma_T \sim 5\text{--}100\text{ mK}$).
3. **Implementation of Simulation Script:** Creation of [`ODE_ring_thermal_noise.m`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/ODE_ring_thermal_noise.m) running single-shot realizations, a $N = 200$ trial Monte Carlo analysis, and a thermal sensitivity parameter sweep.
4. **In-Depth Theoretical Analysis:** Detailed mathematical and physical explanations of the resulting waveforms, frequency response distortions, uncertainty bands, parameter clusters, and the two operating regimes (model-limited vs. thermal-noise-dominated).

---

## 2. Theoretical Framework & Mathematical Equations

### 2.1 Target First-Order ODE
The device implements a first-order Linear Time-Invariant (LTI) optical ODE:
$$\frac{dy(t)}{dt} + a_0 y(t) = \frac{dx(t)}{dt} + b_0 x(t)$$
with spectral transfer function (Coupled-Mode Theory):
$$T_{\text{ODE}}(\omega) = \frac{j\omega + b_0}{j\omega + a_0}$$

### 2.2 Physical MRR Through-Port Transfer Function
The full periodic add-drop MRR model with two couplers is:
$$H_{\text{MRR}}(\omega) = \frac{r_1 - r_2 a_{\text{loss}} \exp\left[-j(\omega \tau_{\text{rt}} + \delta\phi_{\text{rt}})\right]}{1 - r_1 r_2 a_{\text{loss}} \exp\left[-j(\omega \tau_{\text{rt}} + \delta\phi_{\text{rt}})\right]}$$
where:
* $\tau_{\text{rt}} = n_g L / c \approx 2.498\text{ ps}$ (Free Spectral Range: $\text{FSR} \approx 400.24\text{ GHz}$).
* $a_{\text{loss}} = \sqrt{1 - \eta}$ is the round-trip amplitude transmission factor ($\alpha = 8\text{ dB/cm}$).
* $r_{1,2} = \sqrt{1 - \kappa_{1,2}}$ are the field transmission coefficients.

### 2.3 Physical Sources of Thermal Noise
1. **Coupler Phase Shifters (Interferometric Couplers 1 & 2):**
   * Directional coupler base split: $\kappa_0 = 0.0441$.
   * Effective power coupling (Eq. 8, Wu et al.):
     $$\kappa_i = 2\kappa_0(1-\kappa_0)\left[1 + \cos(\Delta\phi_{i,\text{nom}} + \delta\phi_i)\right]$$
   * Thermo-optic phase error:
     $$\delta\phi_i = \frac{2\pi}{\lambda_0} \frac{dn_{\text{eff}}}{dT} L_{bi} \, \delta T_{hi}$$
     where $\frac{dn_{\text{eff}}}{dT} = 1.86 \times 10^{-4}\text{ K}^{-1}$, $L_{b1} = 116.80\ \mu\text{m}$ (Heater 1), and $L_{b2} = 47.12\ \mu\text{m}$ (Heater 2).
   * This perturbs the external decay rates $\gamma_{e1}, \gamma_{e2}$ and drifts ODE coefficients:
     $$a_0 = \gamma_i + \gamma_{e1} + \gamma_{e2}, \quad b_0 = \gamma_i + \gamma_{e2} - \gamma_{e1}$$

2. **Ring Cavity Resonance Detuning & Thermal Cross-Talk:**
   * Heat leakage from microheaters ($\chi = 6\%$) and substrate TEC drift induce temperature variations in the ring:
     $$\delta T_{\text{cavity}} = \delta T_{\text{ambient}} + \chi \left(\frac{\delta T_1 + \delta T_2}{2}\right)$$
   * Round-trip phase error along circumference $L = 178.98\ \mu\text{m}$:
     $$\delta\phi_{\text{rt}} = \frac{2\pi}{\lambda_0} \frac{dn_{\text{eff}}}{dT} L \, \delta T_{\text{cavity}}$$
   * Resonance frequency shift:
     $$\Delta f_{\text{res}} = \frac{\delta\phi_{\text{rt}}}{2\pi \tau_{\text{rt}}} \approx 8.6\text{ GHz/K} \times \delta T_{\text{cavity}}$$

---

## 3. Monte Carlo Methodology

The Monte Carlo routine evaluates statistical sensitivity over $N_{\text{trials}} = 200$ runs:

1. **Stochastic Sampling:** Sample $\delta T_1, \delta T_2 \sim \mathcal{N}(0, \sigma_{T,\text{heater}}^2)$ and $\delta T_{\text{ambient}} \sim \mathcal{N}(0, \sigma_{T,\text{ring}}^2)$.
2. **Phase & Device Mapping:** Compute $\Delta\phi_1^{(k)}, \Delta\phi_2^{(k)}, \kappa_1^{(k)}, \kappa_2^{(k)}$, and $\delta\phi_{\text{rt}}^{(k)}$.
3. **Spectral Propagation:** Evaluate $H_{\text{MRR}}^{(k)}(\omega)$ and filter the input spectrum $X(\omega) = \mathcal{F}\{x(t)\}$ via FFT.
4. **Time-Domain Field & Intensity:** Inverse FFT yields $y^{(k)}(t)$ and power $I^{(k)}(t) = |y^{(k)}(t)|^2$.
5. **Statistical Metrics:**
   * Ensemble mean waveform $\mu_I(t)$ and variance $\sigma_I(t)$.
   * Normalized Mean Square Error (NMSE):
     $$\text{NMSE}^{(k)} = \frac{\sum_t \left|I_{\text{ODE,ideal}}(t) - I^{(k)}(t)\right|^2}{\sum_t \left|I_{\text{ODE,ideal}}(t)\right|^2}$$

---

## 4. Summary of Analysis Figures

| Figure | File | Description & Physical Insight |
| :--- | :--- | :--- |
| **Figure 1** | [`fig1_time_domain_thermal.png`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/fig1_time_domain_thermal.png) | **Time-Domain Output & Error:** Demonstrates split-lobe differentiation $|dx/dt|^2$. Thermal errors cause peak asymmetry and fill the central notch. |
| **Figure 2** | [`fig2_baseband_spectral_thermal.png`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/fig2_baseband_spectral_thermal.png) | **Baseband Spectral Response (Zoomed $\pm 20\text{ GHz}$):** Shows the resonance notch detuning ($\Delta f_{\text{res}} \approx +127.5\text{ MHz}$) and the shifted phase flip. |
| **Figure 3** | [`fig3_monte_carlo_thermal.png`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/fig3_monte_carlo_thermal.png) | **Monte Carlo Statistics ($N=200$):** (a) $\pm 1\sigma$ and $\pm 2\sigma$ pulse confidence bands; (b) NMSE histogram; (c) Anti-correlated $(a_0, b_0)$ parameter drift cluster. |
| **Figure 4** | [`fig4_sensitivity_sweep_thermal.png`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/fig4_sensitivity_sweep_thermal.png) | **Accuracy Sensitivity Sweep ($\sigma_T \in [1, 100]\text{ mK}$):** Transition from model-limited regime ($-60.3\text{ dB}$) to thermal-noise-dominated regime. |

---

## 5. Quantitative Benchmark Comparison

| Metric | Ideal Target ODE | Nominal MRR (Noiseless) | Standard Lab Control ($\sigma_T = 20\text{ mK}$) | Open-Loop / Uncooled ($\sigma_T = 50\text{ mK}$) |
| :--- | :--- | :--- | :--- | :--- |
| **Resonance Detuning $\Delta f_{\text{res}}$** | $0\text{ MHz}$ | $0\text{ MHz}$ | $\pm 125\text{--}250\text{ MHz}$ | $\pm 400\text{--}800\text{ MHz}$ |
| **Coupling $\kappa_1$** | $0.0800$ | $0.0800$ | $0.0800 \pm 0.0003$ | $0.0800 \pm 0.0008$ |
| **Coupling $\kappa_2$** | $0.0400$ | $0.0400$ | $0.0400 \pm 0.0001$ | $0.0400 \pm 0.0004$ |
| **ODE Coefficient $a_0$** | $3.145 \times 10^{10}\text{ rad/s}$ | $3.145 \times 10^{10}\text{ rad/s}$ | $(3.145 \pm 0.005) \times 10^{10}$ | $(3.145 \pm 0.015) \times 10^{10}$ |
| **ODE Coefficient $b_0$** | $-1.919 \times 10^{9}\text{ rad/s}$ | $-1.919 \times 10^{9}\text{ rad/s}$ | $(-1.92 \pm 0.04) \times 10^{9}$ | $(-1.92 \pm 0.11) \times 10^{9}$ |
| **Mean NMSE** | $0$ (Exact) | **$-60.3\text{ dB}$** ($9.4 \times 10^{-7}$) | **$-47.9\text{ dB}$** ($1.6 \times 10^{-5}$) | **$-41.7\text{ dB}$** ($6.8 \times 10^{-5}$) |
| **Accuracy Penalty** | -- | *Baseline* | **$+12.4\text{ dB}$** | **$+18.6\text{ dB}$** |

---

## 6. How to Run the Scripts

To execute the original baseline simulation:
```matlab
ODE_ring
```

To execute the thermal noise Monte Carlo simulation and regenerate all figures:
```matlab
ODE_ring_thermal_noise
```
To adjust thermal noise parameters, edit lines 57–59 of [`ODE_ring_thermal_noise.m`](file:///C:/Users/vitto/OneDrive/Desktop/Uni2025_2026/Photonic_computing/ODE_MRR_solver/ODE_ring_thermal_noise.m#L57-L59):
```matlab
sigma_T_heater = 0.020;  % [K] Heater temperature control accuracy (e.g. 0.005 for high-end TEC)
sigma_T_ring   = 0.015;  % [K] Ring cavity substrate temperature stability
thermal_crosstalk = 0.06; % 6% thermal leakage from microheaters to ring
```
