# Morton CRAE Evapotranspiration Calculator

## Overview
This Python implementation calculates daily evapotranspiration using the Morton Complementary Relationship Areal Evapotranspiration (CRAE) method. The calculator processes meteorological data to estimate potential evapotranspiration (ETP) and actual evapotranspiration (ET) values for long-term periods. The sample data base data is based on 13 years of collected data but can be adapted for any long-term period.

## Data Input Format
The calculator requires daily meteorological data in an Excel workbook with the following columns:

| Column | Variable | Unit | Description |
|--------|----------|------|-------------|
| Date   | t        | YYYY-MM-DD | Date of observation |
| Tmax   | T_max    | °C    | Maximum daily temperature |
| Tavg   | T_avg    | °C    | Average of maximum and minimum temperature |
| TD     | T_d      | °C    | Dew point temperature |
| S      | S        | -     | Sunshine ratio (dimensionless) |
| n      | n_days   | d     | Number of days |
| i      | month    | 1-12  | Month number |

## Key Constants

### Site-Specific Constants
- **Altitude (H)**: 730 m
- **Latitude (φ)**: 55.21664°
- **Long-term mean precipitation (PA)**: 397.1 mm

### CRAE Empirical Constants
- **b1**: 14.0 W/m²
- **b2**: 1.20 W/m²
- **Emissivity (ε)**: 0.92
- **Stefan-Boltzmann constant (σ)**: 5.67×10⁻⁸ W·m⁻²·K⁻⁴
- **b0**: 1.00

## Mathematical Formulations

### 1. Pressure Ratio

$$
\frac{p}{p_s} = \left( \frac{288 - 0.0065H}{288} \right)^{5.256}
$$

Where:  
$p$ = Atmospheric pressure at station (mbar)  
$p_s$ = Sea level pressure (mbar)  
$H$ = Altitude (m)


### 2. Dry-season Snow-free Clear-sky Albedo

$$
a_{zd} = 0.26 - 0.00012 P_A (p/p_s)^{0.5} \left[1 + \frac{\Phi}{42} + \left( \frac{\Phi}{42} \right)^2 \right]
$$

$$
0.11 < a_{zd} < 0.17
$$

Where:  
$a_{zd}$ = Zenith value of dry-season snow-free clear-sky albedo  
$P_A$ = Long-term mean precipitation (mm)  
$p/p_s$ = Pressure ratio (dimensionless)  
$\Phi$ = Latitude (°)


### 3. Saturation Vapour Pressure

$$
e_s(T) = 6.11 \exp\left( \frac{\alpha T}{T + \beta} \right)
$$

Where:  
$e_s$ = Saturation vapour pressure (mbar)  
$T$ = Temperature (°C)  
$\alpha, \beta$ = Empirical constants


### 4. Slope of Saturation Vapour Pressure Curve

$$
\Delta = \frac{\alpha \beta e_s(T)}{(T + \beta)^2}
$$

Where:  
$\Delta$ = Slope of saturation vapour pressure curve (mbar/°C)  
$e_s$ = Saturation vapour pressure (mbar)


### 5. Solar Geometry

#### Solar Declination

$$
\theta = 23.2 \sin\left( \frac{2\pi}{365} (29.5m - 94) \right)
$$

Where:  
$\theta$ = Solar declination (°)  
$m$ = Month number


#### Noon Zenith Distance

$$
\cos(Z) = \cos(\varphi - \theta)
$$

Where:  
$Z$ = Noon zenith distance (°)  
$\varphi$ = Latitude (°)  
$\theta$ = Solar declination (°)


#### Sun Angle

$$
\cos(\omega) = 1 - \frac{\cos(Z)}{\cos(\varphi) - \cos(\theta)}
$$

Where:  
$\omega$ = Sun angle (°)  
$Z$ = Noon zenith distance (°)  
$\varphi$ = Latitude (°)  
$\theta$ = Solar declination (°)


### 6. Extra-Atmospheric Radiation

$$
G_E = \frac{1354}{\eta^2} \left( \frac{\omega}{180} \right) \cos(z)
$$

Where:  
$G_E$ = Extra-atmospheric radiation (W/m²)  
$\eta$ = Sun radius vector  
$\omega$ = Sun angle (°)  
$z$ = Average angular zenith distance (°)


### 7. Clear-sky Albedo Chain

$$
c_0 = \max(0, v - v_D)
$$

$$
a_z = a_{zz} + (1 - c_0^2)(0.34 - a_{zz})
$$

$$
a_0 = a_z \frac{\exp(1.08) - \left[2.16 \cos(Z^\circ) + \sin Z\right] \exp(0.012Z)}{1.473(1 - \sin Z)}
$$

*Note: Cosine of Z in degrees, as used in the Python implementation.*


Where:  
$c_0$ = Cloud cover factor  
$v$ = Vapour pressure (mbar)  
$v_D$ = Dew point vapour pressure (mbar)  
$a_z$ = Zenith value of clear-sky albedo  
$a_0$ = Clear-sky albedo  
$Z$ = Noon zenith distance (°)


### 8. Precipitable Water Vapour

$$
W = \frac{v_D}{0.49 + T/129}
$$

Where:  
$W$ = Precipitable water vapour (mm)  
$v_D$ = Dew point vapour pressure (mbar)  
$T$ = Air temperature (°C)

### 9. Turbidity Coefficient

$$
c_1 = \max(0, \min(21 - T, 5))
$$

$$
j = (0.5 + 2.5 \cos^2 z) \cdot e^{c_1 (p/p_s - 1)}
$$

Where:  
$j$ = Turbidity coefficient  
$T$ = Air temperature (°C)  
$z$ = Average angular zenith distance (°)  
$p/p_s$ = Pressure ratio (dimensionless)


### 10. Transmittance

$$
\tau = \exp\left( -0.089 \left(\frac{p/p_s}{\cos z}\right)^{0.75} - 0.083 \left(\frac{j}{\cos z}\right)^{0.90} - 0.029 \left(\frac{W}{\cos z}\right)^{0.60} \right)
$$

Where:  
$\tau$ = Transmittance  
$p/p_s$ = Pressure ratio  
$z$ = Average angular zenith distance (°)  
$j$ = Turbidity coefficient  
$W$ = Precipitable water vapour (mm)


### 11. Clear-sky Global Radiation

$$
G_0 = G_E \tau \left[1 + \left(1 - \frac{\tau}{\tau_a}\right)(1 + a_0 \tau)\right]
$$

Where:  
$G_0$ = Clear-sky global radiation (W/m²)  
$G_E$ = Extra-atmospheric radiation (W/m²)  
$\tau$ = Transmittance  
$\tau_a$ = Absorption component  
$a_0$ = Clear-sky albedo


### 12. Incident Global Radiation

$$
G = S G_0 + (0.08 + 0.30 S)(1 - S) G_E
$$

Where:  
$G$ = Incident global radiation (W/m²)  
$S$ = Sunshine ratio  
$G_0$ = Clear-sky global radiation (W/m²)  
$G_E$ = Extra-atmospheric radiation (W/m²)


### 13. Mean Albedo

$$
a = a_0 S + a_z (1 - S) \left(1 - \frac{Z}{330} \right)
$$


Where:  
$a$ = Mean albedo  
$a_0$ = Clear-sky albedo  
$S$ = Sunshine ratio  
$a_z$ = Zenith value of clear-sky albedo


### 14. Cloud Radiation Effect

$$
c_2 = \max(0, \min(10(v_D/v - S - 0.42), 1.0))
$$

$$
\rho = 0.18 \left[(1 - c_2)(1 - S)^2 + c_2 (1 - S)^{0.5}\right] \cdot \frac{1}{p/p_s}
$$

Where:  
$\rho$ = Cloud radiation effect  
$v$ = Vapour pressure (mbar)  
$v_D$ = Dew point vapour pressure (mbar)  
$S$ = Sunshine ratio  
$p/p_s$ = Pressure ratio


### 15. Net Long-wave Radiation

$$
B = \max\left(\epsilon \sigma T_K^4 \left[1 - (0.71 + 0.007 v_D \cdot p/p_s)(1 + \rho)\right], 0.05 \epsilon \sigma T_K^4\right)
$$

$$
\text{B}_\text{min} = 0.05 \cdot \varepsilon \cdot \sigma \cdot T_k^4
$$


Where:  
$B$ = Net long-wave radiation (W/m²)  
$\epsilon$ = Emissivity (0.92)  
$\sigma$ = Stefan-Boltzmann constant (5.67×10⁻⁸ W m⁻² K⁻⁴)  
$T_K$ = Air temperature (K = °C + 273.15)  
$v_D$ = Dew point vapour pressure (mbar)  
$p/p_s$ = Pressure ratio  
$\rho$ = Cloud radiation effect


### 16. Stability, Transfer, and Iterative Solution for Potential ET

> **Note:** All equations are written using LaTeX syntax for scientific clarity and precision.

---

### Stability and Transfer Coefficients

Let:

- $T$ = air temperature (°C)  
- $v$, $v_D$ = vapour pressures (mbar)  
- $\Delta$ = slope of saturation vapour pressure curve  
- $R_{Tc}$ = net radiation  
- $p/p_s$ = pressure ratio  

#### Psychrometric Constant

$$
\gamma_p = \frac{\gamma_{ps}}{p/p_s}
$$

where

$$
\gamma_{ps} =
\begin{cases}
0.66, & T \geq 0^\circ\text{C} \\
\frac{0.66}{1.15}, & T < 0^\circ\text{C}
\end{cases}
$$

---

#### Stability Factor

$$
\zeta = \min\left( \frac{1}{\max\left( 0.28 \left( 1 + \frac{v_D}{v} \right ) + \frac{ \Delta R_{Tc} }{ \gamma_p (p/p_s)^{-0.5} b_0 f_z (v - v_D) },\ 1 \times 10^{-6} \right)},\ 1.0 \right)
$$


where

$$
f_z =
\begin{cases}
28.0, & T \geq 0^\circ\text{C} \\
28.0 \times 1.15, & T < 0^\circ\text{C}
\end{cases}
\qquad
b_0 = 1.0
$$

---

#### Vapour Transfer Coefficient

$$
f_T = (p/p_s)^{-0.5} \cdot \frac{f_z}{\zeta}
$$

---

#### Heat Transfer Coefficient

$$
\lambda = \gamma_p + \frac{4 \epsilon \sigma (T + 273.15)^3}{f_T}
$$

---

#### Latent Heat of Vaporization or Sublimation

$$
\Delta H_{\text{vap}} =
\begin{cases}
28.5, & T \geq 0^\circ\text{C} \\
28.5 \times 1.15, & T < 0^\circ\text{C}
\end{cases}
$$

---

### Iterative Solution for Equilibrium Surface Temperature ($T_P$)

Start with initial values:

$$
T'_p = T, \quad v'_p = v, \quad \Delta'_P = \Delta
$$

Update equations:

$$
\delta T_p =
\frac{
\frac{R_T}{f_T} + v_D - v'_p + \lambda (T - T'_p)
}{
\Delta'_P + \lambda
}
$$

$$
T_P = T'_p + \delta T_p
$$

$$
v_P = 6.11 \cdot \exp \left( \frac{\alpha T_P}{T_P + \beta} \right)
$$

$$
\Delta_P = \frac{\alpha \beta v_P}{(T_P + \beta)^2}
$$

Repeat until:

$$
|\delta T_p| \leq 0.01^\circ\text{C}
$$

---

### Final Evapotranspiration Calculations

#### Potential ET at Air Temperature

$$
ET_P = R_T - \lambda \cdot f_T (T_P - T)
$$

#### Net Radiation at Surface Temperature $T_P$

$$
R_{TP} = ET_P + \gamma_p f_T (T_P - T)
$$

#### Wet-Environment Areal ET

$$
ET_W = b_1 + b_2 \left(1 + \frac{\gamma_p}{\Delta_P} \right)^{-1} R_{TP}
$$

#### Actual Areal ET (Complementary Relationship)

$$
ET = 2 ET_W - ET_P
$$

#### Conversion to Millimetres

$$
\text{mm} = \left(\text{W m}^{-2} \right) \cdot \frac{n_{\text{days}}}{\Delta H_{\text{vap}}}
$$

---

### Constants

- $b_1 = 14.0$ W/m²  
- $b_2 = 1.20$ W/m²  
- $n_{\text{days}}$ = number of days in the averaging period

---

### 18. Workflow & Formula Dependency Flowcharts

Below are two conceptual flowcharts showing the dependency of all major formulae in the Morton CRAE calculation workflow. The first flowchart is more specific while the second flowchart is a more generalized overview of the workflow of this code.

```mermaid
graph TD
    A[Input Data: T, TD, S, n, i] --> B[Calculate p/ps, azd, v, vD, Δ]
    B --> C[Solar Geometry: θ, Z, ω, z, η, GE]
    C --> D[Albedos: az, ao]
    D --> E[W, j]
    E --> F[Transmittance: τ, τa]
    F --> G[Clear-sky Global Radiation: G0]
    G --> H[Incident Global Radiation: G]
    H --> I[Mean Albedo: a]
    I --> J[Cloud Radiation Effect: ρ]
    J --> K[Net Longwave Radiation: B]
    K --> L[Net Radiation: RT]
    L --> M[Stability/Transfer: ζ, fT, λ, γp, ΔHvap]
    M --> N[Iterative TP Solution: TP, vP, ΔP]
    N --> O[Final ET Calculation: ETP, ET, RT mm]

```

This is the second flowchart that is more generalized.

```mermaid
flowchart TD
    A["Input Data (Excel)"] --> B["Site Constants"]
    B --> C["Meteorological Calculations"]
    C --> D["Radiation Balance"]
    D --> E["Evapotranspiration Calculations"]
```

---

## Data Sources
The meteorological data used in this calculator was collected from:

1. Environment and Climate Canada
   - Daily temperature and precipitation data
   - Station: Grande Airport (for proxy data)
   - Period: January 2012 to March 2025

2. National Research Council Canada
   - Solar geometry calculations
   - Sun calculator service

## Implementation Details
The implementation uses the following Python packages:
- `pandas`: For data processing and Excel file handling
- `numpy`: For mathematical operations and array handling
- `matplotlib`: For visualization capabilities

## Usage
1. Configure the input file path in the code:
   ```python
   XLS_PATH = "New Input Data for Python Morton CRAE Method.xlsm"
   ```

2. Ensure your Excel file contains the required columns with correct headers

3. Run the script to calculate evapotranspiration values

## References
- Pamula, A., & Szyk, M. (2024). [Magnus-Tetens equation accuracy analysis.](https://www.researchgate.net/publication/385816714_Development_of_an_Optimized_Non-Linear_Model_for_Precise_Dew_Point_Estimation_in_Variable_Environmental_Conditions) Journal of Meteorology.
- Environment and Climate Canada. (2025). [Climate Data Online](https://climate.weather.gc.ca/)
- National Research Council Canada. (2025). [Sun Calculator](https://nrc.canada.ca/en/research-development/products-services/software-applications/sun-calculator/)
- Government of Alberta. (April 2013). [Evaporation and Evapotranspiration in Alberta](https://open.alberta.ca/dataset/46557167-126e-4bb3-b84f-615ead212b3f/resource/93686041-152d-400d-854e-b12d4d3a5481/download/8938.pdf)
