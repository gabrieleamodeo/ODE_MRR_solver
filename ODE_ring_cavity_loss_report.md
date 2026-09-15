# Silicon Photonic ODE Solver: Cavity Propagation Loss Sensitivity & Drop-Port Leakage Analysis

**Reference**: Wu et al., *"Compact tunable silicon photonic differential-equation solver for general linear time-invariant systems"*, Optics Express (2014).  
**Script Implementation**: [`ODE_ring_cavity_loss_analysis.m`](file:///Users/gabrieleamodeo/Documents/polimi/photonic/ODE_MRR_solver/ODE_ring_cavity_loss_analysis.m)

---

## 1. Executive Summary

This report analyzes the physical performance, operational stability, and optical power distribution of an all-optical first-order ordinary differential equation (ODE) solver implemented on a silicon photonic add-drop microring resonator (MRR). Specifically, this investigation models:
1. **Fabrication Sensitivity to Waveguide Loss ($\alpha$)**: Sweeping optical propagation loss across a realistic silicon-on-insulator (SOI) fabrication tolerance window ($\alpha \in [4.0, 14.0]\text{ dB/cm}$) under open-loop (fixed coupler heaters) operation.
2. **Coupling Regimes & Transition to Differentiator Mode**: Demonstrating the analytical transition from the **over-coupled** regime ($b_0 < 0$) to the **critically coupled** regime ($b_0 = 0$, pure optical differentiator) and the **under-coupled** regime ($b_0 > 0$).
3. **Drop-Port Leakage Quantification**: Tracking the secondary drop waveguide, which acts as an optical energy drain and low-pass smoothing channel.
4. **Solver Precision Degradation (NMSE)**: Quantifying the Normalized Mean Square Error between the realized physical through-port optical field and the nominal target ODE mathematical solution.

---

## 2. Device Architecture & Mathematical Foundations

### 2.1 Add-Drop Microring Resonator ODE Solver
The solver is configured around a racetrack microring cavity of circumference $L = 178.98\ \mu\text{m}$ coupled to two bus waveguides:
* **Bus/Through Port**: Serves as the computational port solving the first-order linear constant-coefficient ODE:
  $$\frac{dy(t)}{dt} + a_0 y(t) = \frac{dx(t)}{dt} + b_0 x(t)$$
* **Drop Port**: Functions as an optical leakage/extraction channel that filters cavity energy into the second waveguide.

```
                   Through Port y_through(t)
                     ^
                     | (1 - kappa1)
Input Pulse x(t) --->+====================+ (Bus Waveguide)
                         | kappa1
                     +---+---+
                     |  MRR  |  Circumference L = 178.98 um
                     | Loss  |  Round-trip time tau_rt = 2.498 ps
                     +---+---+
                         | kappa2
                     +===+================+ (Drop Waveguide)
                     |
                     v
                   Drop Port y_drop(t) (Optical Leakage Channel)
```

### 2.2 Coupled-Mode Theory (CMT) Formulation
From temporal Coupled-Mode Theory, the cavity energy amplitude $a(t)$ obeys:
$$\frac{da(t)}{dt} = \left( j(\omega - \omega_0) - (\gamma_i + \gamma_{e1} + \gamma_{e2}) \right) a(t) - j\sqrt{2\gamma_{e1}} s_{\text{in}}(t)$$
where:
* **Optical Carrier Angular Frequency**: $\omega_0 = \frac{2\pi c}{\lambda_0}$ with $\lambda_0 = 1550.391\text{ nm}$.
* **Cavity Quality Factors**:
  $$Q_i = -\frac{\omega_0 n_g L}{c \ln(1 - \eta)}, \quad Q_{e1} = -\frac{\omega_0 n_g L}{c \ln(1 - \kappa_1)}, \quad Q_{e2} = -\frac{\omega_0 n_g L}{c \ln(1 - \kappa_2)}$$
* **Decay Rates** ($\gamma = \frac{\omega_0}{2Q}$):
  $$\gamma_i = \frac{c \, \alpha_{\text{lin}}}{2 n_g}, \quad \gamma_{e1} = -\frac{c \ln(1 - \kappa_1)}{2 n_g L}, \quad \gamma_{e2} = -\frac{c \ln(1 - \kappa_2)}{2 n_g L}$$
* **ODE Coefficients**:
  $$a_0 = \gamma_i + \gamma_{e1} + \gamma_{e2} = \frac{\omega_0}{2 Q_{\text{total}}}$$
  $$b_0 = \gamma_i + \gamma_{e2} - \gamma_{e1}$$

The baseband target ODE transfer function is:
$$T_{\text{ODE}}(\omega) = \frac{j\omega + b_0}{j\omega + a_0}$$

---

## 3. Optical Transfer Functions

### 3.1 Through-Port Transfer Function
Accounting for the exact periodic multi-pass interference inside the microring:
$$H_{\text{through}}(\omega) = \frac{r_1 - r_2 a_{\text{loss}} e^{-j\omega \tau_{\text{rt}}}}{1 - r_1 r_2 a_{\text{loss}} e^{-j\omega \tau_{\text{rt}}}}$$
where:
* $r_1 = \sqrt{1 - \kappa_1}$, $r_2 = \sqrt{1 - \kappa_2}$
* $\tau_{\text{rt}} = \frac{n_g L}{c} \approx 2.498\text{ ps}$ (Free Spectral Range $\text{FSR} \approx 400.24\text{ GHz}$)
* $a_{\text{loss}} = \sqrt{1 - \eta} = \exp\left(-\frac{\alpha_{\text{lin}} L}{2}\right)$

### 3.2 Drop-Port Transfer Function (Leakage Channel)
The optical field coupled into the second waveguide after traversing half the cavity round-trip ($\tau_{\text{rt}}/2$):
$$H_{\text{drop}}(\omega) = \frac{-\sqrt{\kappa_1 \kappa_2 a_{\text{loss}}} e^{-j\omega \tau_{\text{rt}} / 2}}{1 - r_1 r_2 a_{\text{loss}} e^{-j\omega \tau_{\text{rt}}}}$$

Around resonance ($\omega \approx 0$), the drop port exhibits a first-order low-pass (integrating) frequency response:
$$H_{\text{drop}}(\omega) \approx \frac{-C_{\text{drop}}}{j\omega + a_0}$$
meaning the drop port functions as an optical integrator/smoother that siphons optical pulse energy away from the through port.

---

## 4. Analytical Critical Coupling Condition & Regimes

At carrier resonance ($\omega = 0$), destructive interference at the through port is dictated by the numerator $r_1 - r_2 a_{\text{loss}}$:
$$b_0 = 0 \iff \gamma_i = \gamma_{e1} - \gamma_{e2} \iff r_1 = r_2 a_{\text{loss}}$$

Solving for the critical linear propagation loss $\alpha_{\text{crit,lin}}$:
$$\exp\left(-\alpha_{\text{crit,lin}} L\right) = \frac{1 - \kappa_1}{1 - \kappa_2} \implies \alpha_{\text{crit,lin}} = -\frac{1}{L} \ln\left(\frac{1 - \kappa_1}{1 - \kappa_2}\right)$$

Evaluating with $\kappa_1 = 0.08$, $\kappa_2 = 0.04$, and $L = 178.98\ \mu\text{m}$:
$$\alpha_{\text{crit}} = \frac{10}{L \ln(10)} \ln\left(\frac{1 - 0.04}{1 - 0.08}\right) \approx 10.3271\text{ dB/cm}$$

### Operational Regimes:
1. **Over-Coupled Regime ($\alpha < 10.33\text{ dB/cm}$, $b_0 < 0$)**:
   - Coupler 1 decay dominates: $\gamma_{e1} > \gamma_i + \gamma_{e2}$.
   - Nominal design ($\alpha_{\text{nom}} = 8.0\text{ dB/cm}$) operates here with $b_0 = -1.92 \times 10^9\text{ rad/s}$.
   - High through-port transmission and moderate notch depth ($-24.3\text{ dB}$).
2. **Critically Coupled Regime ($\alpha \approx 10.33\text{ dB/cm}$, $b_0 \approx 0$)**:
   - Internal loss exactly balances external coupling: $\gamma_i = \gamma_{e1} - \gamma_{e2}$.
   - The ODE becomes $\frac{dy}{dt} + a_0 y = \frac{dx}{dt}$, which for pulse bandwidths below $a_0$ reduces to a **pure optical differentiator**: $y(t) \propto \frac{dx(t)}{dt}$.
   - Complete destructive cancellation at resonance ($|H_{\text{through}}(0)|^2 \to -\infty\text{ dB}$).
   - The through-port temporal intensity splits into a double-lobed pulse with a zero center null.
3. **Under-Coupled Regime ($\alpha > 10.33\text{ dB/cm}$, $b_0 > 0$)**:
   - Cavity loss dominates: $\gamma_i > \gamma_{e1} - \gamma_{e2}$.
   - The notch depth degrades, transmission increases at DC, and the zero-crossing is lost.

---

## 5. Quantitative Metric Summary Table

The table below summarizes the calculated metrics for the nominal design alongside the three representative operating points investigated in the script:

| Operating Condition | Waveguide Loss $\alpha$ | Intrinsic Rate $\gamma_i$ | Realized $a_0$ | Realized $b_0$ | Coupling Regime | Static Notch Depth | Through NMSE vs Target | Through Energy | Drop Leakage | Cavity Dissipation |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Nominal Reference** | **$8.00\text{ dB/cm}$** | **$6.60 \times 10^9\text{ rad/s}$** | **$3.15 \times 10^{10}\text{ rad/s}$** | **$-1.92 \times 10^9\text{ rad/s}$** | Over-coupled | **$-24.29\text{ dB}$** | **$-69.38\text{ dB}$** | **$29.09\%$** | **$39.22\%$** | **$31.69\%$** |
| **Point A (Under-Loss)** | $5.00\text{ dB/cm}$ | $4.12 \times 10^9\text{ rad/s}$ | $2.90 \times 10^{10}\text{ rad/s}$ | $-4.39 \times 10^9\text{ rad/s}$ | Deep over-coupled | $-16.39\text{ dB}$ | $-22.57\text{ dB}$ | $32.86\%$ | $44.62\%$ | $22.53\%$ |
| **Point B (Critical Cpl.)** | $10.33\text{ dB/cm}$ | $8.52 \times 10^9\text{ rad/s}$ | $3.34 \times 10^{10}\text{ rad/s}$ | $\approx 0\text{ rad/s}$ | Differentiator | $-297.15\text{ dB}$ | $-26.59\text{ dB}$ | $27.10\%$ | $35.68\%$ | $37.22\%$ |
| **Point C (Over-Loss)** | $12.00\text{ dB/cm}$ | $9.90 \times 10^9\text{ rad/s}$ | $3.48 \times 10^{10}\text{ rad/s}$ | $+1.38 \times 10^9\text{ rad/s}$ | Under-coupled | $-28.03\text{ dB}$ | $-22.42\text{ dB}$ | $26.06\%$ | $33.43\%$ | $40.52\%$ |

---

## 6. Generated Visualizations

The script [`ODE_ring_cavity_loss_analysis.m`](file:///Users/gabrieleamodeo/Documents/polimi/photonic/ODE_MRR_solver/ODE_ring_cavity_loss_analysis.m) automatically exports two high-resolution verification figures:

### Figure 1: `fig1_loss_metrics_sweep.png`
Consists of four subplots tracking the loss sweep from $4.0$ to $14.0\text{ dB/cm}$:
* **(a) Realized ODE Coefficients**: Illustrates the monotonic rise of $a_0(\alpha)$ and $b_0(\alpha)$ with loss, clearly identifying the sign-flip of $b_0$ across $\alpha_{\text{crit}} = 10.33\text{ dB/cm}$.
* **(b) Static Extinction Ratio / Notch Depth**: Demonstrates the singular notch plunge (down to theoretical zero) at $\alpha_{\text{crit}}$, where through-port destructive cancellation is complete.
* **(c) Through-Port Solver Error (NMSE)**: Evaluates the sensitivity of the output intensity. While the nominal design matches the ideal target with an intrinsic approximation error of $-69.4\text{ dB}$, uncompensated fabrication drifts of $\pm 2-3\text{ dB/cm}$ degrade NMSE to $\approx -22\text{ dB}$.
* **(d) Optical Pulse Energy Budget Breakdown**: Displays the complementary partition of input optical pulse energy into Through Transmission, Drop Leakage, and Cavity Absorption/Scattering.

### Figure 2: `fig2_loss_time_domain_waveforms.png`
Examines the 10 Gb/s Gaussian pulse ($45\text{ ps}$ FWHM) dynamics across three distinct operating points:
* **Point A ($\alpha = 5.0\text{ dB/cm}$)**: Over-coupled through-port intensity and strong drop leakage ($44.62\%$).
* **Point B ($\alpha = 10.33\text{ dB/cm}$)**: Differentiator mode showing the clean splitting of the input pulse into the characteristic symmetric two-lobe profile with zero intensity at the pulse center.
* **Point C ($\alpha = 12.0\text{ dB/cm}$)**: Under-coupled output where loss increases cavity dissipation to over $40.5\%$.

---

## 7. Key Physical & Engineering Insights

1. **Passive Differentiator Synthesis via Loss Tuning**:
   A standard through-port ODE solver can be directly turned into a high-bandwidth optical differentiator by adjusting the cavity loss to $\alpha_{\text{crit}} = 10.33\text{ dB/cm}$ (or dynamically tuning the coupling ratios $\kappa_1, \kappa_2$ via on-chip microheaters such that $r_1 = r_2 a_{\text{loss}}$).
2. **Drop-Port as a Primary Energy Sink**:
   Across all realistic loss values ($4$ to $14\text{ dB/cm}$), the drop port continually leaks between **$30\%$ and $45\%$ of the total incident optical pulse energy**. This confirms that the drop port is the largest single source of optical insertion loss for the through-port solver.
3. **Necessity of Closed-Loop Coupler Calibration**:
   Fabrication variations in waveguide core geometry or surface roughness directly alter $\alpha$. Because open-loop operation shifts the realized ODE coefficient $b_0$ by several gigaradians per second, active phase shifter tuning on the MZI couplers (as modeled in [`ODE_ring_thermal_noise.m`](file:///Users/gabrieleamodeo/Documents/polimi/photonic/ODE_MRR_solver/ODE_ring_thermal_noise.m)) is essential to re-balance $\kappa_1$ and $\kappa_2$ and preserve target ODE solver accuracy.

---

## 8. Usage Instructions

To execute the simulation in MATLAB:
```matlab
% Navigate to the project directory
cd /Users/gabrieleamodeo/Documents/polimi/photonic/ODE_MRR_solver

% Run the cavity loss analysis
ODE_ring_cavity_loss_analysis
```
The script will run the 101-point sweep, display the formatted metrics in the command window, render Figure 1 and Figure 2, and save them as PNG images in the working directory.
