import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, linear_model, typeA_uncertainty

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 17" / "figures"
REPORT = BASE_DIR / "lab 17" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# =========================
# Constants and parameters
# =========================
R_used = 1974

lambda_red = 635
lambda_yel = 590
lambda_gre = 525
lambda_blu = 460

e_charge = 1.602176634e-19
c_light   = 2.99792458e8
h_true    = 6.62607015e-34

# =========================
# Helper functions
# =========================
def compute_r_squared(x, a, b, y):
    y_predicted = linear_model(x, a, b)
    residuals   = y - y_predicted
    ss_res      = np.sum(residuals**2)
    ss_tot      = np.sum((y - np.mean(y))**2)
    return 1 - (ss_res / ss_tot)

def weighted_linear_fit(x, y):
    popt, pcov = curve_fit(linear_model, x, y)
    slope, intercept = popt
    slope_err, intercept_err = np.sqrt(np.diag(pcov))
    cov_sb    = pcov[0, 1]
    r_squared = compute_r_squared(x, slope, intercept, y)
    return slope, intercept, slope_err, intercept_err, r_squared, cov_sb

def compute_relative_error(measured, true):
    return (measured - true) / true * 100.0

def find_tangent_region(V, I, window=5):
    dIdV  = np.gradient(I, V)
    idx   = np.argmax(dIdV)
    start = max(0, idx - window)
    end   = min(len(V), idx + window)
    return start, end

# =========================
# Load data
# =========================
df = pd.read_excel(BASE_DIR / "lab 17" / "data.xlsx")
df["I"] = df["I"] * 0.001

# =========================
# Subexperiment A: Resistor
# =========================
df_A = df[df["subexp_id"] == "A"]
V_A  = df_A["V"].values
I_A  = df_A["I"].values

s1, b1, s1_err, b1_err, r1, _ = weighted_linear_fit(V_A, I_A)

R         = 1.0 / s1
R_err     = s1_err / (s1**2)
rel_err_R = compute_relative_error(R, R_used)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(V_A, I_A, label="Data Point")
V_fit = np.linspace(min(V_A), max(V_A), 100)
if b1 < 0:
    plt.plot(V_fit, linear_model(V_fit, s1, b1), color='black', linestyle='--',
             label=f'Best-fit line: y = {s1:.4g} x - {abs(b1):.4g} (R² = {r1:.4f})')
else:
    plt.plot(V_fit, linear_model(V_fit, s1, b1), color='black', linestyle='--',
             label=f'Best-fit line: y = {s1:.4g} x + {b1:.4g} (R² = {r1:.4f})')
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title("Subexperiment for Resistor: I-V Curve")
plt.legend()
plt.savefig(FIGURES_DIR / "A_IV_fit.png", dpi=300)
plt.close()

# =========================
# Subexperiment B: Diode
# =========================
df_B = df[df["subexp_id"] == "B"]

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.plot(0.001 * df_B["V"], df_B["I"], "o-", label="Data Point")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.title("Subexperiment for Silicon Diode: I-V Curve")
plt.legend()
plt.savefig(FIGURES_DIR / "B_IV.png", dpi=300)
plt.close()

# =========================
# Subexperiment C: LED
# =========================
df_C = df[df["subexp_id"] == "C"]

lambda_map = {
    "red":    lambda_red,
    "yellow": lambda_yel,
    "green":  lambda_gre,
    "blue":   lambda_blu,
}

color_caption = {
    "blue":   "藍色LED在不同壓降下流經的電流",
    "red":    "紅色LED在不同壓降下流經的電流",
    "yellow": "黃色LED在不同壓降下流經的電流",
    "green":  "綠色LED在不同壓降下流經的電流",
}

VF_list      = []
VF_err_list  = []
lambda_list  = []
color_list   = []
h_list       = []
h_err_list   = []
c_data_tables = []

for color in df_C["color"].unique():
    df_color = df_C[df_C["color"] == color]
    V = df_color["V"].values
    I = df_color["I"].values

    idx_sort = np.argsort(V)
    V = V[idx_sort]
    I = I[idx_sort]

    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.axhline(0, color="gray", alpha=0.7, linewidth=0.8)
    plt.plot(V, I, "o-", label="Data Point")
    plt.xlabel("Voltage (V)")
    plt.ylabel("Current (A)")
    plt.title(f"Subexperiment for {color} LED: I-V Curve")

    start, end = find_tangent_region(V, I)
    V_fit_c = V[start:end]
    I_fit_c = I[start:end]

    s_tan, b_tan, s_tan_err, b_tan_err, r_tan, cov_sb = weighted_linear_fit(V_fit_c, I_fit_c)

    VF     = -b_tan / s_tan
    VF_err = np.sqrt(
        (b_tan_err / s_tan)**2 +
        (b_tan * s_tan_err / s_tan**2)**2 -
        2 * (b_tan / s_tan**3) * cov_sb
    )

    K         = e_charge * 10**(-9) * lambda_map[color] / c_light
    h         = VF * K
    h_err     = VF_err * K
    rel_err_h = compute_relative_error(h, h_true)

    V_line = np.linspace(VF, max(V), 100)
    if b_tan < 0:
        plt.plot(V_line, linear_model(V_line, s_tan, b_tan), color='black', linestyle='--',
                 label=f'Tangent-fit line: y = {s_tan:.4g} x - {abs(b_tan):.4g} (R² = {r_tan:.4f})')
    else:
        plt.plot(V_line, linear_model(V_line, s_tan, b_tan), color='black', linestyle='--',
                 label=f'Tangent-fit line: y = {s_tan:.4g} x + {b_tan:.4g} (R² = {r_tan:.4f})')
    plt.legend()
    plt.savefig(FIGURES_DIR / f"C_IV_{color}.png", dpi=300)
    plt.close()

    VF_list.append(VF)
    VF_err_list.append(VF_err)
    lambda_list.append(lambda_map[color])
    color_list.append(color)
    h_list.append(pm(h, h_err))
    h_err_list.append(percent(rel_err_h))

    c_data_tables.append((
        pd.DataFrame({
            "Voltage (\\unit{\\V})":  df_color["V"].map(lambda x: f"${x:.4g}$"),
            "Current (\\unit{\\mA})": df_color["I"].map(lambda x: f"${1000*x:.4g}$"),
        }),
        color_caption.get(color, f"{color.capitalize()}色LED在不同壓降下流經的電流"),
    ))

# =========================
# VF vs 1/lambda (all class)
# =========================
VF_mean_list   = []
VF_unc_list    = []
lambda_vf_list = []
VF_data_list   = []
h_result_list  = []
h_error_list   = []

df_vf = pd.read_excel(BASE_DIR / "lab 17" / "turnonvoltage.xlsx")

def remove_outliers(data):
    data  = np.array(data, dtype=float)
    Q1    = np.percentile(data, 25)
    Q3    = np.percentile(data, 75)
    IQR   = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    return data[(data >= lower) & (data <= upper)]

for lam, group in df_vf.groupby("lambda"):
    VF_values = group["VF"].to_numpy()
    VF_data_list.append(VF_values)

    VF_clean = remove_outliers(VF_values)
    mean, u  = typeA_uncertainty(VF_clean)

    lambda_vf_list.append(lam)
    VF_mean_list.append(mean)
    VF_unc_list.append(u)

    K     = e_charge * 10**(-9) * lam / c_light
    h     = mean * K
    h_err = u * K
    h_result_list.append(pm(h, h_err))
    h_error_list.append(percent(compute_relative_error(h, h_true)))

fig, ax = plt.subplots(figsize=(8, 5))
ax.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
ax.boxplot(VF_data_list, labels=[f"{lam}" for lam in lambda_vf_list], vert=False)
ax.set_ylabel(r"$\lambda$ (nm)")
ax.set_xlabel(r"$V_F$ (V)")
fig.savefig(FIGURES_DIR / "boxplot.png", dpi=300)
plt.close()

lambda_arr  = np.array(lambda_vf_list, dtype=float)
VF_mean_arr = np.array(VF_mean_list)
VF_unc_arr  = np.array(VF_unc_list)
inv_lambda  = 10**9 / lambda_arr

popt, pcov = curve_fit(
    linear_model,
    inv_lambda,
    VF_mean_arr,
    sigma=VF_unc_arr,
    absolute_sigma=True,
)
s2, b2         = popt
s2_err, b2_err = np.sqrt(np.diag(pcov))
r2             = compute_r_squared(inv_lambda, s2, b2, VF_mean_arr)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(inv_lambda, VF_mean_arr, yerr=VF_unc_arr, fmt='o', capsize=3, label="Data Point")
x_fit = np.linspace(min(inv_lambda), max(inv_lambda), 100)
if b2 < 0:
    plt.plot(x_fit, linear_model(x_fit, s2, b2), color='black', linestyle='--',
             label=f'Best-fit line: y = {s2:.4g} x - {abs(b2):.4g} (R² = {r2:.4f})')
else:
    plt.plot(x_fit, linear_model(x_fit, s2, b2), color='black', linestyle='--',
             label=f'Best-fit line: y = {s2:.4g} x + {b2:.4g} (R² = {r2:.4f})')
plt.xlabel(r"$1/\lambda$ ($\mathrm{m}^{-1}$)")
plt.ylabel(r"$V_F$ (V)")
plt.legend()
plt.savefig(FIGURES_DIR / "VF_vs_inv_lambda.png", dpi=300)
plt.close()

K_planck       = e_charge / c_light
h_planck       = s2 * K_planck
h_err_planck   = s2_err * K_planck
rel_err_planck = compute_relative_error(h_planck, h_true)

# =========================
# Output tables (in main.tex order)
# =========================

df_to_latex(
    pd.DataFrame({
        "Voltage (\\unit{\\V})":  df_A["V"].map(lambda x: f"${x:.4g}$"),
        "Current (\\unit{\\mA})": df_A["I"].map(lambda x: f"${1000*x:.4g}$"),
    }),
    caption="電阻在不同壓降下流經的電流",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "slope $s$",
            "intercept",
            "$R_{\\text{obs}}=1/s$ (\\unit{\\ohm})",
            "$R_{\\text{multi}}$ (\\unit{\\ohm})",
            "Error",
        ],
        "Value": [
            pm(s1, s1_err),
            pm(b1, b1_err),
            pm(R, R_err),
            f"${R_used}$",
            percent(rel_err_R),
        ],
    }),
    caption="電阻值的測量結果",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Voltage (\\unit{\\mV})":  df_B["V"].map(lambda x: f"${x:.4g}$"),
        "Current (\\unit{\\mA})":  df_B["I"].map(lambda x: f"${1000*x:.4g}$"),
    }),
    caption="矽二極體在不同壓降下流經的電流，負號代表逆向偏壓",
    path=REPORT,
)

for df_table, caption in c_data_tables:
    df_to_latex(df_table, caption=caption, path=REPORT)

df_to_latex(
    pd.DataFrame({
        "Color": [c.capitalize() for c in color_list],
        "Turn-on Voltage $V_F$ (\\unit{\\V})": [
            pm(v, u) for v, u in zip(VF_list, VF_err_list)
        ],
        "$\\lambda$ (\\unit{\\nm})": [f"${lam}$" for lam in lambda_list],
        "$h_{\\text{obs}}=eV_F\\lambda/c$ (\\unit{\\J\\s})": h_list,
        "$h_{\\text{ref}}$ (\\unit{\\J\\s})": ["$6.626\\times 10^{-34}$"] * len(color_list),
        "Error": h_err_list,
    }),
    caption="導通電壓和普朗克常數的計算結果",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "$\\lambda$ (\\unit{\\nm})":         [f"${lam}$" for lam in lambda_vf_list],
        "$V_F$ (\\unit{\\V})":               [pm(v, u) for v, u in zip(VF_mean_list, VF_unc_list)],
        "$h_{\\text{obs}}$ (\\unit{\\J\\s})": h_result_list,
        "Error":                             h_error_list,
    }),
    caption="去除異常值後，使用各個數據集計算普朗克常數的結果",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "slope",
            "intercept",
            "$h_{\\text{obs}}$ (\\unit{\\J\\s})",
            "$h_{\\text{ref}}$ (\\unit{\\J\\s})",
            "Error",
        ],
        "Value": [
            pm(s2, s2_err),
            pm(b2, b2_err),
            pm(h_planck, h_err_planck),
            "$6.626\\times 10^{-34}$",
            percent(rel_err_planck),
        ],
    }),
    caption="由上圖回歸直線斜率計算普朗克常數的結果",
    path=REPORT,
)
