# Morton CRAE Daily Evapotranspiration Estimator
# Author: Martin Edwini-Bonsu
# Language: Python 3.13.2

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# 0.  ▌USER SETTINGS▐
# ─────────────────────────────────────────────────────────────────────────────
XLS_PATH   = "New Input Data for Python Morton CRAE Method.xlsm"   # full path to your Excel workbook
DATA_SHEET = "daily"                    # sheet containing daily observations
# Expected columns in DATA_SHEET:
#  Date | Tmax | T | TD | S | n | i
#       |  °C  |°C| °C | – | d | 1-12
#
#  If your file uses different headers, edit COLS below.
# Tmax = Maximum temperature
# Tavg = Average temperature
# TD = Dew point temperature
# S = Sunshine ratio
# n = Number of days
# i = Month (1-12)
COLS = dict(date="Date", tmax="Tmax", tavg="tavg", tdew="TD",
            sun_ratio="S", ndays="n", month="i")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  ▌SITE CONSTANTS (modify for other sites)▐
# ─────────────────────────────────────────────────────────────────────────────
H   = 730              # altitude, m
phi = 55.21664         # latitude, °
PA  = 397.1            # long-term mean precip, mm

# CRAE empirical constants
b1, b2   = 14.0, 1.20
epsilon  = 0.92
sigma    = 5.67e-8     # W m-2 K-4
b0       = 1.00

# ─────────────────────────────────────────────────────────────────────────────
# 2.  ▌HELPER FUNCTIONS▐
#    Each follows the numbering of your spec for easy audit.
# ─────────────────────────────────────────────────────────────────────────────
def p_ps(alt_m):
    """Eq 2 – pressure ratio p/ps""" ## Compute ratio of atmospheric pressure at the station to that at sea level with the pressure correction equation 
    # for the standard atmosphere
    return ((288 - 0.0065*alt_m)/288) ** 5.256

# Estimate the zenith value of the dry-season snow-free clear-sky albedo
def azd(pa_mm, p_over_ps, lat_deg):
    """Eq 3 – dry-season, snow-free clear-sky albedo"""
    azd_raw = 0.26 - 0.00012*pa_mm*(p_over_ps**0.5)*(1 + lat_deg/42 + (lat_deg/42)**2)
    return np.clip(azd_raw, 0.11, 0.17)

# Compute the saturation vapour pressure at dew point temperature and the saturation vapour pressure
def sat_vap_press(t_c, alpha, beta):
    """Eq 3/4 – saturation vapour pressure (mbar)"""
    return 6.11*np.exp(alpha*t_c/(t_c+beta))

# Computer the slope of the saturation vapour pressure curve
def slope_svp(t_c, v_mbar, alpha, beta):
    """Eq 5 – slope of the SVP curve Δ (mbar °C-1)"""
    return (alpha*beta*v_mbar)/(t_c+beta)**2

# Compute the various angles and functions leading up to an estimate of the extra-atmospheric global radiation
# theta is the declination of the sun (in degrees)
# omega(w) is the number of degrees the earth rotates between sunrise and noon
# Z and z are the noon and average angular zenith distances of the sun, respectively
# eta (η) is the radius vector of the sun
def solar_geometry(month, lat, p_over_ps, v, vD):
    """Steps 6-7 – returns dict with θ, Z, ω, z, η, GE, azz, az, ao"""
    theta = 23.2*np.sin(np.deg2rad(29.5*month - 94))       # Eq 6
    cosZ  = np.cos(np.deg2rad(lat - theta))                # Eq 7
    cosZ  = np.where(cosZ < 0.001, 0.001, cosZ)            # Eq 7a
    cosω  = 1 - cosZ / (np.cos(np.deg2rad(lat)) - np.cos(np.deg2rad(theta)))  # Eq 8
    cosω  = np.clip(cosω, -1.0, 1.0)  # keep strictly in domain for arccos
    ω     = np.rad2deg(np.arccos(cosω))

    # Avoid /0 if ω==0 (polar night) by small floor:
    # Estimate the zenith value of snow-free clear-sky albedo (azz), the zenith value of clear-sky albedo (az), and the clear-sky albedo (ao)
    ω     = np.where(ω == 0, 1e-6, ω)
    term  = (180/np.pi)*(np.sin(np.deg2rad(ω))/ω) - 1
    cosz  = cosZ + term*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(theta))   # Eq 9
    η     = 1 + (1/60)*np.sin(np.deg2rad(29.5*month - 106))                # Eq 10
    GE    = (1354/η**2)*(ω/180)*cosz                                       # Eq 11
    # Clear-sky albedo chain
    azz   = None  # will be set outside because depends on azd
    return dict(theta=theta, Z=np.rad2deg(np.arccos(cosZ)),
                ω=ω, z=np.rad2deg(np.arccos(cosz)), η=η, GE=GE, azz=azz)

def albedos(azz, v, vD, Z):
    """Eqs 12-15 – returns az, ao"""
    c0  = v - vD
    c0  = np.clip(c0, 0, 1)                                # Eq 13a
    az  = azz + (1 - c0**2)*(0.34 - azz)                   # Eq 14
    # Eq 15 rewritten safely (Z in degrees)
    num = np.exp(1.08) - (2.16*np.cos(np.deg2rad(Z)) + np.sin(np.deg2rad(Z)))*np.exp(0.012*Z)
    den = 1.473*(1 - np.sin(np.deg2rad(Z)))
    ao  = az*num/den
    return az, ao

# Estimate precipitable water vapour W in millimetres and a turbidity coefficient
def precipitable_W(vD, T):
    """Eq 16 – precipitable water vapour W (mm)"""
    return vD/(0.49 + T/129)

def turbidity_j(T, z_deg, p_ps_ratio):
    """Eq 18 – turbidity coefficient j"""
    c1 = np.clip(21 - T, 0, 5)                             # Eq 17 & 17a
    return (0.5 + 2.5*np.cos(np.deg2rad(z_deg))**2) * np.exp(c1*(p_ps_ratio - 1))

# Compute the transmittancy of clear skies to direct beam solar radiation from an equation formulated by Brooks
def transmittance_tau(p_ps_ratio, z_deg, j, W):
    """Eq 9 – τ"""
    cosz = np.cos(np.deg2rad(z_deg))
    term1 = -0.089*(p_ps_ratio/cosz)**0.75
    term2 = -0.083*(j/cosz)**0.90
    term3 = -0.029*(W/cosz)**0.60
    return np.exp(term1 + term2 + term3)

# Estimate the part of transmittancy of clear skies to direct beam solar radiation that is the result of absorption (tau_a)

def tau_a(j, W, z_deg):
    """Eq 10 – τa (absorption component)"""
    cosz = np.cos(np.deg2rad(z_deg))
    return np.exp(-0.0415*(j/cosz)**0.90 - (0.0029)**0.50*(W/cosz)**0.30)

#Compute the clear-sky global radiation Go (W m-2)
def clear_sky_G0(GE, τ, τa, ao):
    """Eq 11-1 – clear-sky global radiation Go (W m-2)"""
    return GE*τ*(1 + (1 - τ/τa)*(1 + ao*τ))

# Compute the incident global radiation (G)
def incident_G(S, Go, GE):
    """Eq 11-2 – incident global radiation G (W m-2)"""
    return S*Go + (0.08 + 0.30*S)*(1 - S)*GE

# EStimate the average albedo (a)
def mean_albedo(a0, S, Z):
    """Eq 12 – average albedo a"""
    return a0*(S + (1 - S)*(1 - Z/330))

# Estimate the proportional increase in atmospheric radiation due to clouds (rho)
def rho_cloud(T, vD, v, S, p_ps_ratio):
    """Eq 13 – ρ, ∝ increase in LW due to clouds"""
    c2 = np.clip(10*(vD/v - S - 0.42), 0, 1.0)
    return 0.18*((1 - c2)*(1 - S)**2 + c2*(1 - S)**0.5) * (1/p_ps_ratio)

# Calculate the net long-wave radiation loss for soil--plant surfaces at air temperature
def B_longwave(T, vD, p_ps_ratio, rho):
    """Eq 14 – net LW radiation loss B (W m-2)"""
    T_k = T + 273.15
    raw = epsilon*sigma*T_k**4*(1 - (0.71 + 0.007*vD*p_ps_ratio)*(1 + rho))
    min_b = 0.05*epsilon*sigma*T_k**4
    return np.maximum(raw, min_b)

# Estimate the net radiation for soil-plant surfaces at air temperature, the stability factor, the vapour transfer coefficient, and the heat transfer coefficient.
def RT_net(G, a, B):
    """Eq 27 – net radiation at air temperature RT"""
    return (1 - a)*G - B

# Additional functions for ζ, fT, λ, iterative TP, etc.
def stability_and_transfer(v, vD, Δ, RTc, T, T_cond, p_ps_ratio):
    """Eqs 29-31 – ζ, fT, λ; handles conditional T≥0 logic inside."""
    # γps and fz switch with temperature sign
    if T_cond >= 0:
        fz, γps = 28.0, 0.66
        ΔHvap   = 28.5        # W-days kg-1
    else:
        fz, γps = 28.0*1.15, 0.66/1.15
        ΔHvap   = 28.5*1.15

    γp   = γps * (1/p_ps_ratio)
    invζ = 0.28*(1 + vD/v) + Δ*RTc/(γp*(p_ps_ratio**-0.5)*b0*fz*(v - vD))
    ζ    = np.minimum(1/np.maximum(invζ, 1e-6), 1.0)   # guard
    fT   = (p_ps_ratio**-0.5)*fz/ζ
    λ    = γp + 4*epsilon*sigma*(T+273.15)**3/fT
    return ζ, fT, λ, γp, ΔHvap

# Choose initial values of T'_p, v'_p, and Δ'_P equal to T_p, v, and Δ 
# and estimate the final values from the following quickly converging iterative 
# solution of the vapour transfer and energy-balance equations:

# Equation (32):
# [δT_p] = [R_n/f_T + v_D - v'_p + λ(T - T'_p)] / (Δ'_p + λ)

# Equation (33):
# T_P = T'_p + [δT_p]

# Equation (34):
# v_p = 6.11 * exp[(α * T_P) / (T_P + β)]

# Equation (35):
# Δ_P = α * β * v_p / (T_P + β)^2

# Equations 32 to 35 are repeated, setting T'_p, v'_p, and Δ'_p equal to the values 
# of T_P, v_P, and Δ_P derived from the preceding iteration until |δT_p| <= 0.01°C. 
# The purpose is to estimate the potential evapotranspiration equilibrium temperature (T_P) 
# from a solution of the vapour transfer and energy-balance equations for a small moist surface.

def iterate_TP(T, vD, RT, fT, λ, *, alpha, beta, tol=0.01, max_iter=50):
    """
    Step 16 – solves for equilibrium surface temperature TP.
    Accepts either a float or a NumPy array for T.
    Returns TP, vP, ΔP  (same shape as T).
    """
    # Ensure TP is a mutable NumPy array
    TP = np.asarray(T, dtype=float).copy()

    for _ in range(max_iter):
        vP  = sat_vap_press(TP, alpha, beta)
        ΔP  = slope_svp(TP, vP, alpha, beta)
        δTP = (RT / fT + vD - vP + λ * (T - TP)) / (ΔP + λ)
        TP += δTP
        if np.all(np.abs(δTP) <= tol):
            break

    vP = sat_vap_press(TP, alpha, beta)
    ΔP = slope_svp(TP, vP, alpha, beta)
    return TP, vP, ΔP

# Estimate the potential evapotranspiration, the net radiation for soil--plant surfaces at the equilibrium temperature and the wet environment areal evapotranspiration
# in which the constants b1 and b2 are 14 W/m^2 and 1.20 W/m^2, respectively.

# Convert the net radiation for soil–plant surfaces at air temperature (R_T), 
# the potential evapotranspiration (E_TP), and the areal evapotranspiration (E_T) 
# from the power units of W m^-2 to the evaporation units of millimetres of depth. 
# This is done by dividing by the latent heat of vaporization or sublimation 
# and multiplying by the number of days.

# The latent heat of vaporization (for T >= 0°C) is:
#   28.5 W-days per kilogram

# The latent heat of sublimation (for T < 0°C) is:
#   28.5 × 1.15 = 32.775 W-days per kilogram

# --- function definition -------------------------------------------
def final_ET(RT, fT, TP, T, γp, ΔP, ΔHvap, ndays, λ):
    """Steps 17-19 – returns ETP, ET, RT (all in mm)."""
    ETP = RT - λ * (fT * (TP - T))
    RTP = ETP + γp * fT * (TP - T)
    ETW = b1 + b2 * (1 + γp / ΔP) ** -1 * RTP
    ET  = 2 * ETW - ETP                     # complementary relationship
    to_mm = ndays / ΔHvap
    return ETP * to_mm, ET * to_mm, RT * to_mm

# ─────────────────────────────────────────────────────────────────────────────
# 3.  ▌LOAD DATA▐
# ─────────────────────────────────────────────────────────────────────────────
COLS = {
    "Date": "date", "Tmax": "tmax", "tavg": "tavg", "TD": "tdew",
    "S": "sun_ratio", "n": "ndays", "i": "month",
}

try:
    raw = pd.read_excel(XLS_PATH, sheet_name=DATA_SHEET)
except FileNotFoundError:
    # ─ demo fallback using one synthetic record ─
    raw = pd.DataFrame({
        "Date": ["2025-07-01"],
        "Tmax": [25.4],
        "tavg": [18.0],
        "TD":   [12.0],
        "S":    [0.62],
        "n":    [1],
        "i":    [7],
    })

# Rename to internal conventions & basic checks
df = raw.rename(columns=COLS)
missing = {"tavg", "tdew", "sun_ratio"} - set(df.columns)
assert not missing, f"Missing columns in sheet: {missing}"

df["date"] = pd.to_datetime(df["date"])

# Static site factors
pps = p_ps(H)
azz_const = azd(PA, pps, phi)

# ─────────────────────────────────────────────────────────────────────────────
# 4.  ▌VECTORISED CALCULATIONS▐
# ─────────────────────────────────────────────────────────────────────────────
# Alpha/Beta arrays depend on sign of T
alpha = np.where(df.tavg >= 0, 17.27, 21.88)
beta  = np.where(df.tavg >= 0, 237.3, 265.5)

# vD, v, Δ
df["vD"] = sat_vap_press(df.tdew, 17.27, 237.3)
df["v"]  = sat_vap_press(df.tavg, alpha, beta)
df["Δ"]  = slope_svp(df.tavg, df.v, alpha, beta)

# Solar geometry
geo = df.apply(lambda r: solar_geometry(r.month, phi, pps, r.v, r.vD), axis=1, result_type="expand")
df = pd.concat([df, geo], axis=1)
df["azz"] = azz_const

# Albedos az, ao
az_ao = df.apply(lambda r: albedos(r.azz, r.v, r.vD, r.Z), axis=1, result_type="expand")
df[["az","ao"]] = az_ao

# W, j
df["W"] = precipitable_W(df.vD, df.tavg)
df["j"] = turbidity_j(df.tavg, df.z, pps)

# τ, τa
df["τ"]  = transmittance_tau(pps, df.z, df.j, df.W)
df["τa"] = tau_a(df.j, df.W, df.z)

# G0, G
df["Go"] = clear_sky_G0(df.GE, df["τ"], df["τa"], df.ao)
df["G"]  = incident_G(df.sun_ratio, df.Go, df.GE)

# a, ρ, B, RT
df["a"]  = mean_albedo(df.ao, df.sun_ratio, df.Z)
df["ρ"]  = rho_cloud(df.tavg, df.vD, df.v, df.sun_ratio, pps)
df["B"]  = B_longwave(df.tavg, df.vD, pps, df["ρ"])
df["RT"] = RT_net(df.G, df.a, df.B)

# ζ, fT, λ, γp, ΔHvap
stab = df.apply(
    lambda r: stability_and_transfer(r.v, r.vD, r["Δ"], r.RT, r.tavg, r.tavg, pps),
    axis=1, result_type="expand")
df[["ζ","fT","λ","γp","ΔHvap"]] = stab

# TP iteration
TP_out = df.apply(
    lambda r: iterate_TP(r.tavg, r.vD, r.RT, r.fT, r.λ,
                         alpha= 17.27 if r.tavg >= 0 else 21.88,
                         beta = 237.3 if r.tavg >= 0 else 265.5),
    axis=1, result_type="expand")
df[["TP","vP","ΔP"]] = TP_out

# Final ETP, ET, RT (mm)
# --- call from df.apply --------------------------------------------
etrip = df.apply(
    lambda r: final_ET(
        r.RT, r.fT, r.TP, r.tavg, r.γp, r.ΔP, r.ΔHvap, r.ndays, r.λ
    ),
    axis=1,
    result_type="expand"
)
df[["ETP_mm", "ET_mm", "RT_mm"]] = etrip.round(3)

# ── NEW: force ET to zero when negative ──
df["ET_mm"] = df["ET_mm"].clip(lower=0)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  ▌RESULTS▐
# ─────────────────────────────────────────────────────────────────────────────
print(df[["date","ETP_mm","ET_mm","RT_mm"]])
df.to_excel('Morton_CRAE_evapotranspiration_estimates (with no negative ET values).xlsx', index=False)

# Optional quick plot
df.plot(x="date", y="ET_mm", marker="o")
plt.ylabel("Areal ET (mm d⁻¹)")
plt.title("Morton CRAE Areal Evapotranspiration")
plt.tight_layout()
plt.show()
