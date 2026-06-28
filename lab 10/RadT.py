import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, linear_model, compute_r_squared

BASE_DIR    = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 10" / "figures"
REPORT      = BASE_DIR / "lab 10" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# Hardcoded data: filament temperature T (K) and sensor output Rad (mV)
T_data = np.array([576.6602192, 825.5902513, 1000.4360, 1140.390916, 1255.576082,
                    1357.080292, 1454.94173, 1523.862293, 1606.3231, 1694.6583,
                    1737.2424, 1796.0901])
Rad_data = np.array([0.4, 1.3, 3.1, 5.6, 8.8, 12.1, 16.4, 20.0, 25.1, 30.1, 35.3, 40.3])

# ==================================================
# Table 5 / T.png: Rad - T, power-law fit Rad = A * T^B + C
# ==================================================
def power_model(T, A, B, C):
    return A * T**B + C

popt, pcov = curve_fit(power_model, T_data, Rad_data, p0=[0.0000001, 4, 0])
A_opt, B_opt, C_opt = popt
u_A, u_B, u_C = np.sqrt(np.diag(pcov))

T_fit = np.linspace(min(T_data), max(T_data), 200)
Rad_fit = power_model(T_fit, *popt)

if C_opt >= 0:
    curve_label = f'Best-fit curve: y = {A_opt:.4g} x^{B_opt:.4g} + {C_opt:.4g}'
else:
    curve_label = f'Best-fit curve: y = {A_opt:.4g} x^{B_opt:.4g} - {abs(C_opt):.4g}'

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(T_data, Rad_data, label="Data Point")
plt.plot(T_fit, Rad_fit, color='black', linestyle='--', label=curve_label)
plt.xlabel("T (K)")
plt.ylabel("Rad (mV)")
plt.title("Rad - T")
plt.legend()
plt.savefig(FIGURES_DIR / "T.png", dpi=300)
plt.close()

df_to_latex(
    pd.DataFrame({
        "Parameter": ["$A$", "$B$", "$C$"],
        "Value": [pm(A_opt, u_A), pm(B_opt, u_B), pm(C_opt, u_C)],
    }),
    caption="上圖中回歸直線的相關資訊",
    path=REPORT,
)

# ==================================================
# Table 6 / T4.png: Rad - T^4, linear fit
# ==================================================
T4_data = T_data**4

popt4, pcov4 = curve_fit(linear_model, T4_data, Rad_data)
slope4, intercept4 = popt4
u_slope4, u_intercept4 = np.sqrt(np.diag(pcov4))
r_squared4 = compute_r_squared(T4_data, slope4, intercept4, Rad_data)

x_fit4 = np.linspace(min(T4_data), max(T4_data), 200)
y_fit4 = linear_model(x_fit4, slope4, intercept4)

if intercept4 >= 0:
    label4 = f'Best-fit line: y = {slope4:.4g} x + {intercept4:.4g} (R² = {r_squared4:.4f})'
else:
    label4 = f'Best-fit line: y = {slope4:.4g} x - {abs(intercept4):.4g} (R² = {r_squared4:.4f})'

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(T4_data, Rad_data, label="Data Point")
plt.plot(x_fit4, y_fit4, color='black', linestyle='--', label=label4)
plt.xlabel("$T^{4}$ (K$^{4}$)")
plt.ylabel("Rad (mV)")
plt.title("Rad - $T^{4}$")
plt.legend()
plt.savefig(FIGURES_DIR / "T4.png", dpi=300)
plt.close()

df_to_latex(
    pd.DataFrame({
        "Parameter": ["slope", "y-intercept"],
        "Value": [pm(slope4, u_slope4), pm(intercept4, u_intercept4)],
    }),
    caption="上圖中回歸直線的相關資訊",
    path=REPORT,
)

# ==================================================
# Table 7 / log.png: log Rad - log T, linear fit
# Theoretical slope = 4 per Stefan-Boltzmann law
# ==================================================
logT   = np.log10(T_data)
logRad = np.log10(Rad_data)

popt_log, pcov_log = curve_fit(linear_model, logT, logRad)
slope_log, intercept_log = popt_log
u_slope_log, u_intercept_log = np.sqrt(np.diag(pcov_log))
r_squared_log = compute_r_squared(logT, slope_log, intercept_log, logRad)

err_slope_log = (slope_log - 4) / 4 * 100

x_fit_log = np.linspace(min(logT), max(logT), 200)
y_fit_log = linear_model(x_fit_log, slope_log, intercept_log)

if intercept_log >= 0:
    label_log = f'Best-fit line: y = {slope_log:.4g} x + {intercept_log:.4g} (R² = {r_squared_log:.4f})'
else:
    label_log = f'Best-fit line: y = {slope_log:.4g} x - {abs(intercept_log):.4g} (R² = {r_squared_log:.4f})'

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(logT, logRad, label="Data Point")
plt.plot(x_fit_log, y_fit_log, color='black', linestyle='--', label=label_log)
plt.xlabel(r"$\log T$")
plt.ylabel(r"$\log$ Rad")
plt.title(r"$\log$ Rad - $\log T$")
plt.legend()
plt.savefig(FIGURES_DIR / "log.png", dpi=300)
plt.close()

df_to_latex(
    pd.DataFrame({
        "Parameter": ["slope", "Error", "y-intercept"],
        "Value": [
            pm(slope_log, u_slope_log),
            percent(err_slope_log),
            pm(intercept_log, u_intercept_log),
        ],
    }),
    caption="上圖中回歸直線的相關資訊，根據史蒂芬-波茲曼定律描述的是整個波段，它應該是一條斜率為四的斜直線",
    path=REPORT,
)
