import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty

BASE_DIR    = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 19" / "figures"
REPORT      = BASE_DIR / "lab 19" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# ==================================================
# Constants
# ==================================================
tau   = 1.005e-3   # RC time constant (s)
C_cap = 1e-9       # capacitance (F)

# Damped RLC (tables 5-6): R = 100+50+48.1 Ω, L = 22 mH
R_damp         = 100 + 50 + 48.1
L_damp         = 22e-3
beta_theo_damp  = R_damp / (2 * L_damp)
omega1_theo_damp = np.sqrt(1 / (C_cap * L_damp) - beta_theo_damp**2)

# Resonance RLC (table 8): R = 100+50+25.8 Ω, L = 10 mH
R_reson          = 100 + 50 + 25.8
L_reson          = 10e-3
Vin_reson        = 1.02   # (Vpp)
beta_theo_reson  = R_reson / (2 * L_reson)
omega1_theo_reson = np.sqrt(1 / (C_cap * L_reson) - 2 * beta_theo_reson**2)
f_reson_theo     = omega1_theo_reson / (2 * np.pi)
q_reson          = R_reson * np.sqrt(C_cap / L_reson)
Q_theo           = omega1_theo_reson / (2 * beta_theo_reson)

# ==================================================
# Load data
# ==================================================
df1 = pd.read_excel(BASE_DIR / "lab 19" / "data1.xlsx")

df2    = pd.read_excel(BASE_DIR / "lab 19" / "data2.xlsx")
f2     = df2["f"].to_numpy(dtype=float)
Vout2  = df2["Vout"].to_numpy(dtype=float)
phase2 = df2["phase"].to_numpy(dtype=float)
idx2   = np.argsort(f2)
f2, Vout2, phase2 = f2[idx2], Vout2[idx2], phase2[idx2]

# Resonant frequency = frequency at peak amplitude
f_reson_obs = float(f2[np.argmax(Vout2)])
A_max       = float(np.max(Vout2))

# Hardcoded damping data
t_peak_us = np.array([15.2, 45.0, 75.6, 106.0, 136.0, 166.0, 197.0])
V_peak    = np.array([4.40, 3.76, 3.20, 2.72, 2.32, 2.00, 1.68])
t_peak_s  = t_peak_us * 1e-6

# ==================================================
# Tables 2, 3, 4 and Gain plots (subexps A, B, C)
# ==================================================
table_captions_gain = {
    "A": "使用正弦波，測量電容壓降，驗證其低通濾波器的屬性",
    "B": "使用正弦波，測量電阻壓降，驗證其高通濾波器的屬性",
    "C": "使用方波，測量電容壓降",
}

for subexp in ["A", "B", "C"]:
    sub   = df1[df1["subexp_id"] == subexp]
    T     = sub["T"].to_numpy(dtype=float)
    Vin_g = sub["Vin"].to_numpy(dtype=float)
    Vout_g = sub["Vout"].to_numpy(dtype=float)
    f_g   = 1.0 / (T * tau)
    Gain  = 20 * np.log10(Vout_g / Vin_g)
    idx_g = np.argsort(f_g)
    f_g, Gain = f_g[idx_g], Gain[idx_g]

    df_to_latex(
        pd.DataFrame({
            "Frequency (\\unit{\\Hz})": [f"${v:.2f}$" for v in f_g],
            "Gain (\\unit{\\decibel})": [f"${v:.2f}$" for v in Gain],
        }),
        caption=table_captions_gain[subexp],
        path=REPORT,
    )

    logf_arr = np.linspace(math.log(min(f_g)), math.log(max(f_g)), 200)
    fs = np.exp(logf_arr)
    if subexp == "A":
        G_theory = 20 * np.log10(
            (1 / (2*np.pi*fs)) / np.sqrt(tau**2 + (1 / (2*np.pi*fs))**2)
        )
    elif subexp == "B":
        G_theory = 20 * np.log10(1 / np.sqrt(1 + (1 / (2*np.pi*tau*fs))**2))
    else:
        G_theory = 20 * np.log10(
            (1 - np.exp(-1 / (fs * 2 * tau))) / (1 + np.exp(-1 / (fs * 2 * tau)))
        )

    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.semilogx(f_g, Gain, 'o', label="Data Point")
    plt.semilogx(fs, G_theory, color='black', linestyle='--', label="Theory prediction")
    plt.xlabel("f (Hz)")
    plt.ylabel("Gain (dB)")
    plt.title("Gain vs f")
    plt.grid(True, which="both")
    plt.legend()
    plt.savefig(FIGURES_DIR / f"Sub-experiment_{subexp}.png", dpi=300)
    plt.close()

# ==================================================
# Table 5: Raw peak data (from hardcoded fitting data)
# ==================================================
df_to_latex(
    pd.DataFrame({
        "Time (\\unit{\\us})":   [f"${v:.3g}$" for v in t_peak_us],
        "Voltage (\\unit{\\V})": [f"${v:.2f}$" for v in V_peak],
    }),
    caption="電容壓降曲線的高峰頂點的時間和電壓",
    path=REPORT,
)

# ==================================================
# Table 6: β and ω₁ results  (was manually created)
# ==================================================
def exp_model(t, A, beta):
    return A * np.exp(-beta * t)

popt, pcov = curve_fit(exp_model, t_peak_s, V_peak)
A_fit, beta_obs = popt
_, beta_err = np.sqrt(np.diag(pcov))

# Quasi period from successive peak time differences
diffs_us            = np.diff(t_peak_us)
period_mean_us, period_u_us = typeA_uncertainty(diffs_us)
omega1_obs = 2 * np.pi / (period_mean_us * 1e-6)
omega1_u   = omega1_obs * (period_u_us / period_mean_us)

err_beta   = (beta_obs   - beta_theo_damp)  / beta_theo_damp  * 100
err_omega1 = (omega1_obs - omega1_theo_damp) / omega1_theo_damp * 100

df_to_latex(
    pd.DataFrame({
        "Parameters": [
            "$\\beta_{\\text{obs}}$ (\\unit{\\per\\s})",
            "$\\beta_{\\text{theo}}$ (\\unit{\\per\\s})",
            "Error",
            "兩相臨高點的時間差 (Quasi) Period (\\unit{\\us})",
            "$\\omega_{1,\\text{obs}}$ (\\unit{\\radian\\per\\s})",
            "$\\omega_{1,\\text{theo}}$ (\\unit{\\radian\\per\\s})",
            "Error",
        ],
        "Values": [
            pm(beta_obs, beta_err),
            f"${beta_theo_damp:.2f}$",
            percent(err_beta),
            pm(period_mean_us, period_u_us),
            pm(omega1_obs, omega1_u),
            f"${omega1_theo_damp:.0f}$",
            percent(err_omega1),
        ],
    }),
    caption="$\\beta$和$\\omega_{1}$的測量結果",
    path=REPORT,
)

# fitting.png
t_sample = np.linspace(min(t_peak_s), max(t_peak_s), 200)
plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(t_peak_s, V_peak, label="Data Point")
plt.plot(t_sample, exp_model(t_sample, A_fit, beta_obs), color='black', linestyle='--',
         label=f"Best-fit curve: y = {A_fit:.2f} exp($-${beta_obs:.0f} x)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.title("Exponential Fitting")
plt.legend()
plt.savefig(FIGURES_DIR / "fitting.png", dpi=300)
plt.close()

# ==================================================
# Table 8: Resonance parameters  (was manually created)
# ==================================================
threshold = A_max / np.sqrt(2)

above     = f2 > f_reson_obs
below     = f2 < f_reson_obs
f_plus    = float(f2[above][np.argmin(np.abs(Vout2[above] - threshold))])
f_minus   = float(f2[below][np.argmin(np.abs(Vout2[below] - threshold))])
delta_f   = f_plus - f_minus
Q_obs     = f_reson_obs / delta_f

err_f_reson = (f_reson_obs - f_reson_theo) / f_reson_theo * 100
err_Q       = (Q_obs - Q_theo) / Q_theo * 100

df_to_latex(
    pd.DataFrame({
        "Parameters": [
            "$\\beta_{\\text{theo}}$ (\\unit{\\per\\s})",
            "$f_{\\text{reson,obs}}$ (\\unit{\\Hz})",
            "$f_{\\text{reson,theo}}$ (\\unit{\\Hz})",
            "Error",
            "$f_{+}$ s.t. $A=A_{\\max}/\\sqrt{2}$ (\\unit{\\Hz})",
            "$f_{-}$ (\\unit{\\Hz})",
            "$\\Delta f$ (\\unit{\\Hz})",
            "$Q_{\\text{obs}}$",
            "$Q_{\\text{theo}}$",
            "Error",
        ],
        "Values": [
            f"${beta_theo_reson:.0f}$",
            f"${f_reson_obs:.0f}$",
            f"${f_reson_theo:.2f}$",
            percent(err_f_reson),
            f"${f_plus:.0f}$",
            f"${f_minus:.0f}$",
            f"${delta_f:.0f}$",
            f"${Q_obs:.2f}$",
            f"${Q_theo:.2f}$",
            percent(err_Q),
        ],
    }),
    caption=(
        f"共振頻率$f_{{\\text{{reson}}}}$和品質因子的測量結果。"
        f"$A_{{\\max}}/\\sqrt{{2}}={A_max:.1f}/\\sqrt{{2}}\\approx{threshold:.2f}$ Vpp"
    ),
    path=REPORT,
)

# ==================================================
# Table 9: Resonance raw data (data2.xlsx)
# ==================================================
df_to_latex(
    pd.DataFrame({
        "Frequency (\\unit{\\kHz})": [f"${v/1000:.3g}$" for v in f2],
        "Ampl (Vpp)":               [f"${v:.3g}$"       for v in Vout2],
        "Phase (\\unit{\\degree})": [f"${v:.3g}$"       for v in phase2],
    }),
    caption="第四個子實驗的數據",
    path=REPORT,
)

# ==================================================
# Resonance plots
# ==================================================

# --- Semilog: all data ---
plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2, 3), useMathText=True)
plt.semilogx(f2, Vout2, 'o-', label="Data Point")
plt.xlabel("f (Hz)")
plt.ylabel("Ampl (Vpp)")
plt.title("Ampl vs f")
plt.grid(True, which="both")
plt.legend()
plt.savefig(FIGURES_DIR / "Vlogf.png", dpi=300)
plt.close()

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2, 3), useMathText=True)
plt.semilogx(f2, phase2, 'o-', label="Data Point")
plt.xlabel("f (Hz)")
plt.ylabel("Phase (degree)")
plt.title("Phase vs f")
plt.grid(True, which="both")
plt.legend()
plt.savefig(FIGURES_DIR / "phaselogf.png", dpi=300)
plt.close()


def make_theory_curves(W_arr):
    Vt = 1 / np.sqrt(q_reson**2 * W_arr**2 + (1 - W_arr**2)**2)
    pt = np.arctan(q_reson / (W_arr - 1 / W_arr))
    pt[pt > 0] -= np.pi
    return Vt, np.rad2deg(pt)


# --- Normalized linear: all data (Vf.png, phasef.png) ---
W2       = f2 / f_reson_obs
W_sample = np.linspace(min(W2), max(W2), 5000)
Vt_all, pt_all = make_theory_curves(W_sample)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.plot(W2, Vout2 / Vin_reson, 'o', label="Data Point")
plt.plot(W_sample, Vt_all, color='black', linestyle='--', label="Theory prediction")
plt.xlabel(r"$f/f_{\mathrm{reson}}$")
plt.ylabel(r"Ampl/$\varepsilon_{0}$")
plt.title(r"Ampl/$\varepsilon_{0}$ vs $f/f_{\mathrm{reson}}$")
plt.grid(True)
plt.legend()
plt.savefig(FIGURES_DIR / "Vf.png", dpi=300)
plt.close()

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.plot(W2, phase2, 'o', label="Data Point")
plt.plot(W_sample, pt_all, color='black', linestyle='--', label="Theory prediction")
plt.xlabel(r"$f/f_{\mathrm{reson}}$")
plt.ylabel("Phase (degree)")
plt.title(r"Phase vs $f/f_{\mathrm{reson}}$")
plt.grid(True)
plt.legend()
plt.savefig(FIGURES_DIR / "phasef.png", dpi=300)
plt.close()

# --- Normalized linear: highest frequency removed (Vfscale.png, phasefscale.png) ---
mask        = f2 < np.max(f2)
W2_zoom     = f2[mask] / f_reson_obs
Vout2_zoom  = Vout2[mask]
phase2_zoom = phase2[mask]
W_zoom      = np.linspace(min(W2_zoom), max(W2_zoom), 5000)
Vt_zoom, pt_zoom = make_theory_curves(W_zoom)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.plot(W2_zoom, Vout2_zoom / Vin_reson, 'o', label="Data Point")
plt.plot(W_zoom, Vt_zoom, color='black', linestyle='--', label="Theory prediction")
plt.xlabel(r"$f/f_{\mathrm{reson}}$")
plt.ylabel(r"Ampl/$\varepsilon_{0}$")
plt.title(r"Ampl/$\varepsilon_{0}$ vs $f/f_{\mathrm{reson}}$ (highest freq removed)")
plt.grid(True)
plt.legend()
plt.savefig(FIGURES_DIR / "Vfscale.png", dpi=300)
plt.close()

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.plot(W2_zoom, phase2_zoom, 'o', label="Data Point")
plt.plot(W_zoom, pt_zoom, color='black', linestyle='--', label="Theory prediction")
plt.xlabel(r"$f/f_{\mathrm{reson}}$")
plt.ylabel("Phase (degree)")
plt.title(r"Phase vs $f/f_{\mathrm{reson}}$ (highest freq removed)")
plt.grid(True)
plt.legend()
plt.savefig(FIGURES_DIR / "phasefscale.png", dpi=300)
plt.close()
