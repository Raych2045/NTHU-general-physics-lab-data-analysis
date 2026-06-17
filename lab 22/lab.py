import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty, linear_model, compute_r_squared

# ==================================================
# CONSTANTS
# ==================================================
lambda_ref = 657e-9
t = 6e-3
n_glass_ref = 1.46
h = 3e-2
p_atm = 101.325

# ==================================================
# CONFIG
# ==================================================
BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 22" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

REPORT = BASE_DIR / "lab 22" / "report.tex"
REPORT.write_text("", encoding="utf-8")

df = pd.read_excel(BASE_DIR / "lab 22" / "data.xlsx")

# ==================================================
# COMMON CALCULATION
# ==================================================
df["delta_read"] = df["final"] - df["initial"]

# ==================================================
# SUBEXPERIMENT A
# ==================================================
df_A = df[df["subexpid"] == "A"].copy()

delta_read_A = df_A["delta_read"].to_numpy(dtype=float)
deltaN_A = df_A["deltaN"].to_numpy(dtype=float)

lambda_A = (2 * delta_read_A / 1000) / (25 * deltaN_A)
lambda_unc = (0.01/1000)/np.sqrt(12)*(2/25/100)
lambda_percent_error = (lambda_A[-1] - lambda_ref) / lambda_ref * 100

vfunc = np.vectorize(lambda x: f"${x:.4g}$")

table_lambda_data = pd.DataFrame({
    "測微頭讀數差 (\\unit{\\mm})": df_A["delta_read"].map(lambda x: f"${x:.3f}$"),
    "N": df_A["deltaN"].map(lambda x: f"${x:.4g}$"),
    "$\\lambda$ (\\unit{\\nm})": vfunc(lambda_A * 1e9),
})
df_to_latex(table_lambda_data, path=REPORT,
            caption="不同條紋數下的雷射波長的計算值")

table_A = pd.DataFrame({
    "Parameter": [
        "$\\lambda_{\\text{obs}}$ (\\unit{\\nm})",
        "$\\lambda_{\\text{ref}}$ (\\unit{\\nm})",
        "Error",
    ],
    "Value": [
        pm(lambda_A[-1] * 1e9, lambda_unc * 1e9),
        f"${lambda_ref * 1e9:.4g}$",
        percent(lambda_percent_error),
    ],
})
df_to_latex(table_A, path=REPORT,
            caption="$N=100$對應的$\\lambda$與參考值的相對誤差，有考慮B類不確定度")

# ==================================================
# SUBEXPERIMENT B
# ==================================================
df_B = df[df["subexpid"] == "B"].copy()

delta_read_B = df_B["delta_read"].to_numpy(dtype=float)
deltaN_B = df_B["deltaN"].to_numpy(dtype=float)

delta_read_B = np.deg2rad(delta_read_B)

aux_B = deltaN_B * lambda_ref
A = 2 * t - aux_B

numerator = A * (1 - np.cos(delta_read_B))
denominator = 2 * t * (1 - np.cos(delta_read_B)) - aux_B

n_glass = numerator / denominator
n_glass_avg, n_glass_unc = typeA_uncertainty(n_glass)

theta_resolution_deg = 1
sigma_theta = np.deg2rad(theta_resolution_deg / np.sqrt(12))

dn_dtheta = (
    A * np.sin(delta_read_B) * denominator
    - numerator * (2 * t * np.sin(delta_read_B))
) / (denominator**2)

sigma_n_theta = np.abs(dn_dtheta) * sigma_theta
sigma_total = np.sqrt(np.sum(sigma_n_theta**2)) / len(sigma_n_theta)

n_glass_percent_error = (n_glass_avg - n_glass_ref) / n_glass_ref * 100

table_glass_data = pd.DataFrame({
    "角度讀數差 (\\unit{\\degree})": df_B["delta_read"].map(lambda x: f"${x:.4g}$"),
    "N": df_B["deltaN"].map(lambda x: f"${x:.4g}$"),
    "$n_{\\text{glass}}$": vfunc(n_glass),
})
df_to_latex(table_glass_data, path=REPORT,
            caption="在此視每組數據皆是獨立的，最後一行四捨五入至小數點下第三位")

table_B = pd.DataFrame({
    "Parameter": [
        "$n_{\\text{glass,obs}}$",
        "$n_{\\text{glass,ref}}$",
        "Error",
    ],
    "Value": [
        pm(n_glass_avg, n_glass_unc),
        f"${n_glass_ref:.4g}$",
        percent(n_glass_percent_error),
    ],
})
df_to_latex(table_B, path=REPORT,
            caption="玻璃折射率測量值與參考值的相對誤差，僅考慮A類不確定度")

# ==================================================
# SUBEXPERIMENT C
# ==================================================
df_C = df[df["subexpid"] == "C"].copy()

delta_read_C = df_C["delta_read"].to_numpy(dtype=float)
deltaN_C = df_C["deltaN"].to_numpy(dtype=float)

delta_n = (lambda_ref * deltaN_C) / (2 * h)

popt, pcov = curve_fit(linear_model, delta_read_C, delta_n)
s, b = popt

x_fit = np.linspace(np.min(delta_read_C), np.max(delta_read_C), 500)
y_fit = linear_model(x_fit, s, b)

n_air = 1 + s * p_atm

sigma_s = np.sqrt(pcov[0, 0])
sigma_b = np.sqrt(pcov[1, 1])
sigma_n_air = p_atm * sigma_s

r_squared = compute_r_squared(delta_read_C, s, b, delta_n)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(delta_read_C, delta_n, label="Data Point")
if b >= 0:
    fit_label = f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {r_squared:.4f})'
else:
    fit_label = f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {r_squared:.4f})'
plt.plot(x_fit, y_fit, color='black', linestyle='--', label=fit_label)
plt.xlabel(r"$\Delta P$ (kPa)")
plt.ylabel(r"$\Delta n$")
plt.title(r"$\Delta n$ vs $\Delta P$")
plt.legend()
plt.savefig(FIGURES_DIR / "delta_n-delta_p.png", dpi=300)
plt.close()

table_air_data = pd.DataFrame({
    "N": df_C["deltaN"].map(lambda x: f"${x:.4g}$"),
    "$\\Delta n$": vfunc(delta_n),
    "$\\Delta P$ (\\unit{\\kPa})": df_C["delta_read"].map(lambda x: f"${x:.1f}$"),
})
df_to_latex(table_air_data, path=REPORT,
            caption="不同壓力差下的空氣折射率差，中間一行只取三位有效數字")

table_C = pd.DataFrame({
    "Parameter": [
        "Slope $s$",
        "Intercept",
        "$n_{\\text{air,obs}}=1+s\\cdot\\qty{101.325}{\\kPa}$",
        "$n_{\\text{air,ref}}$",
    ],
    "Value": [
        pm(s, sigma_s),
        pm(b, sigma_b),
        pm(n_air, sigma_n_air),
        "$1.00025{-}1.00035$",
    ],
})
df_to_latex(table_C, path=REPORT, caption="外推回歸直線而計算出的空氣折射率")
