import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty, linear_model, compute_r_squared

# =========================
# constants
# =========================
g = 9.8         # m/s^2, standard gravity
M_ref = 209.2e-3  # kg, rotating-body mass (balance reading), fixed in both experiments
r1 = 0.1200       # m, rotation radius fixed in Experiment 1 (also one of Experiment 2's radii)
m_ref = 55.34e-3  # kg, hanging-mass balance reading, fixed in Experiment 2

M_ref_g = M_ref * 1000
r1_cm = r1 * 100
m_ref_g = m_ref * 1000

T_TYPEB_S = 0.010 / np.sqrt(12)  # s, Type B uncertainty of T from the 10 ms timer resolution

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "Lab3" / "figures"
REPORT = BASE_DIR / "Lab3" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# =========================
# load data
# =========================
df1_raw = pd.read_excel(BASE_DIR / "Lab3" / "data.xlsx", sheet_name="Exp1")
df2_raw = pd.read_excel(BASE_DIR / "Lab3" / "data.xlsx", sheet_name="Exp2")


def period_stats(raw, group_col):
    rows = []
    for key, sub in raw.groupby(group_col):
        T_mean_s, u_A_s = typeA_uncertainty(sub["T_s"].to_numpy(dtype=float))
        u_T_s = np.sqrt(u_A_s**2 + T_TYPEB_S**2)
        rows.append({group_col: key, "T_mean_s": T_mean_s, "u_T_s": u_T_s})
    return pd.DataFrame(rows).sort_values(group_col).reset_index(drop=True)


def raw_table(raw, group_col, header):
    keys = sorted(raw[group_col].unique())
    return pd.DataFrame({
        header: [f"${k:.2f}$" for k in keys],
        "Raw periods $T$ (\\unit{\\s})": [
            ", ".join(f"{t:.2f}" for t in raw.loc[raw[group_col] == k, "T_s"])
            for k in keys
        ],
    })


stats1 = period_stats(df1_raw, "m_g")
stats2 = period_stats(df2_raw, "r_cm")

# =========================
# Experiment 1: hanging mass m -> period T -> inferred mass
# =========================
stats1["m_kg"] = stats1["m_g"] / 1000
stats1["m_calc_kg"] = (4 * np.pi**2 * M_ref * r1 / g) / stats1["T_mean_s"]**2
stats1["u_m_calc_kg"] = (4 * np.pi**2 * M_ref * r1 / g) * (2 * stats1["u_T_s"] / stats1["T_mean_s"]**3)
stats1["error_pct"] = (stats1["m_calc_kg"] - stats1["m_kg"]) / stats1["m_kg"] * 100

table1 = pd.DataFrame({
    "Mass $m$ (\\unit{\\g})": stats1["m_g"].map(lambda x: f"${x:.2f}$"),
    "Period $T$ (\\unit{\\s})": [pm(t, u) for t, u in zip(stats1["T_mean_s"], stats1["u_T_s"])],
    "Inferred mass $m_{\\text{obs}}$ (\\unit{\\g})": [pm(m * 1000, u * 1000) for m, u in zip(stats1["m_calc_kg"], stats1["u_m_calc_kg"])],
    "Error": stats1["error_pct"].map(percent),
})
df_to_latex(table1, path=REPORT,
            caption=f"Inferred hanging mass from the measured period, with $M=\\qty{{{M_ref_g:.1f}}}{{\\g}}$ and $r=\\qty{{{r1_cm:.2f}}}{{\\cm}}$ fixed")

x1 = stats1["m_kg"].to_numpy()
y1 = 1 / stats1["T_mean_s"].to_numpy()**2
sigma_y1 = 2 * stats1["u_T_s"].to_numpy() / stats1["T_mean_s"].to_numpy()**3

popt1, pcov1 = curve_fit(linear_model, x1, y1, sigma=sigma_y1, absolute_sigma=True)
s1, b1 = popt1
u_s1, u_b1 = np.sqrt(np.diag(pcov1))
r_squared1 = compute_r_squared(x1, s1, b1, y1)

x1_fit = np.linspace(min(x1), max(x1), 100)
y1_fit = linear_model(x1_fit, *popt1)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(x1, y1, yerr=sigma_y1, fmt="o", capsize=3, label="Data Point")
if b1 >= 0:
    fit_label1 = f'Best-fit line: y = {s1:.4g} x + {b1:.4g} (R² = {r_squared1:.4f})'
else:
    fit_label1 = f'Best-fit line: y = {s1:.4g} x - {abs(b1):.4g} (R² = {r_squared1:.4f})'
plt.plot(x1_fit, y1_fit, color='black', linestyle='--', label=fit_label1)
plt.xlabel(r"Hanging mass $m$ (kg)")
plt.ylabel(r"$T^{-2}$ (s$^{-2}$)")
plt.legend()
plt.savefig(FIGURES_DIR / "T-2-m.png", dpi=300)
plt.close()

M_calc_kg = g / (4 * np.pi**2 * r1 * s1)
u_M_calc_kg = M_calc_kg * abs(u_s1 / s1)
error_M_pct = (M_calc_kg - M_ref) / M_ref * 100

table2 = pd.DataFrame({
    "Quantity": [
        "Slope $s$ (\\unit{\\per\\kg\\per\\s\\squared})",
        "$M_{\\text{obs}}$ (\\unit{\\g})",
        "$M_{\\text{ref}}$ (\\unit{\\g})",
        "Error",
    ],
    "Value": [
        pm(s1, u_s1),
        pm(M_calc_kg * 1000, u_M_calc_kg * 1000),
        f"${M_ref_g:.1f}$",
        percent(error_M_pct),
    ],
})
df_to_latex(table2, path=REPORT,
            caption="Linear regression of $T^{-2}$ against $m$ in Experiment 1, and the resulting estimate of $M$")

# =========================
# Experiment 2: rotation radius r -> period T
# =========================
stats2["r_m"] = stats2["r_cm"] / 100

table3 = pd.DataFrame({
    "Rotation radius $r$ (\\unit{\\cm})": stats2["r_cm"].map(lambda x: f"${x:.2f}$"),
    "Period $T$ (\\unit{\\s})": [pm(t, u) for t, u in zip(stats2["T_mean_s"], stats2["u_T_s"])],
})
df_to_latex(table3, path=REPORT,
            caption=f"Measured period at different rotation radii, with $m$ fixed at an unknown value and $M=\\qty{{{M_ref_g:.1f}}}{{\\g}}$")

x2 = stats2["r_m"].to_numpy()
y2 = stats2["T_mean_s"].to_numpy()**2
sigma_y2 = 2 * stats2["T_mean_s"].to_numpy() * stats2["u_T_s"].to_numpy()

popt2, pcov2 = curve_fit(linear_model, x2, y2, sigma=sigma_y2, absolute_sigma=True)
s2, b2 = popt2
u_s2, u_b2 = np.sqrt(np.diag(pcov2))
r_squared2 = compute_r_squared(x2, s2, b2, y2)

x2_fit = np.linspace(min(x2), max(x2), 100)
y2_fit = linear_model(x2_fit, *popt2)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(x2, y2, yerr=sigma_y2, fmt="o", capsize=3, label="Data Point")
if b2 >= 0:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x + {b2:.4g} (R² = {r_squared2:.4f})'
else:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x - {abs(b2):.4g} (R² = {r_squared2:.4f})'
plt.plot(x2_fit, y2_fit, color='black', linestyle='--', label=fit_label2)
plt.xlabel(r"Rotation radius $r$ (m)")
plt.ylabel(r"$T^{2}$ (s$^{2}$)")
plt.legend()
plt.savefig(FIGURES_DIR / "T2-r.png", dpi=300)
plt.close()

m_calc2_kg = 4 * np.pi**2 * M_ref / (g * s2)
u_m_calc2_kg = m_calc2_kg * abs(u_s2 / s2)
error_m2_pct = (m_calc2_kg - m_ref) / m_ref * 100

table4 = pd.DataFrame({
    "Quantity": [
        "Slope $s$ (\\unit{\\s\\squared\\per\\m})",
        "$m_{\\text{obs}}$ (\\unit{\\g})",
        "$m_{\\text{ref}}$ (\\unit{\\g})",
        "Error",
    ],
    "Value": [
        pm(s2, u_s2),
        pm(m_calc2_kg * 1000, u_m_calc2_kg * 1000),
        f"${m_ref_g:.2f}$",
        percent(error_m2_pct),
    ],
})
df_to_latex(table4, path=REPORT,
            caption="Linear regression of $T^{2}$ against $r$ in Experiment 2, and the resulting estimate of $m$")

# =========================
# Cross-check: predict m from Experiment 1's T^-2-m fit, using Experiment 2's
# period at the radius matching Experiment 1's fixed r
# =========================
row_cross = stats2[np.isclose(stats2["r_cm"], r1_cm)].iloc[0]
T_cross_s = row_cross["T_mean_s"]
u_T_cross_s = row_cross["u_T_s"]

y_cross = 1 / T_cross_s**2
u_y_cross = 2 * u_T_cross_s / T_cross_s**3

m_cross_kg = (y_cross - b1) / s1
u_m_cross_kg = np.sqrt(
    (u_y_cross / s1)**2
    + (m_cross_kg**2 * u_s1**2 + u_b1**2) / s1**2
    + 2 * m_cross_kg * pcov1[0, 1] / s1**2
)
error_cross_pct = (m_cross_kg - m_ref) / m_ref * 100

table_cross = pd.DataFrame({
    "Quantity": [
        "$m_{\\text{obs}}$ (\\unit{\\g})",
        "$m_{\\text{ref}}$ (\\unit{\\g})",
        "Error",
    ],
    "Value": [
        pm(m_cross_kg * 1000, u_m_cross_kg * 1000),
        f"${m_ref_g:.2f}$",
        percent(error_cross_pct),
    ],
})
df_to_latex(table_cross, path=REPORT,
            caption=f"Mass predicted from Experiment 1's fit line, using Experiment 2's period at $r=\\qty{{{r1_cm:.2f}}}{{\\cm}}$")

# =========================
# Raw data (appendix)
# =========================
df_to_latex(raw_table(df1_raw, "m_g", "Mass $m$ (\\unit{\\g})"), path=REPORT,
            caption="Raw period measurements for each hanging mass in Experiment 1")
df_to_latex(raw_table(df2_raw, "r_cm", "Rotation radius $r$ (\\unit{\\cm})"), path=REPORT,
            caption="Raw period measurements for each rotation radius in Experiment 2")
