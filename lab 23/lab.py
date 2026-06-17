import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, linear_model

# =========================
# constants
# =========================
d = 10000 / 6  # nm
lambda_ref_3 = 656.3  # nm
lambda_ref_4 = 486.1  # nm
lambda_ref_5 = 434    # nm
lambda_ref_6 = 410.2  # nm
alpha = np.pi / 3     # rad

h_ref = 6.62607015e-34  # J s
c = 2.99792458e8        # m/s
e = 1.602176634e-19     # C
m_e = 9.1093837e-31     # kg
eps0 = 8.854187817e-12  # F/m
m_p = 1.67262192e-27    # kg

mu = 1 / (1/m_e + 1/m_p)
R_ref = mu * e**4 / (8 * eps0**2 * h_ref**3 * c)

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 23" / "figures"
REPORT = BASE_DIR / "lab 23" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# =========================
# load data
# =========================
df = pd.read_excel(BASE_DIR / "lab 23" / "data.xlsx")

# =========================
# angle processing
# =========================
df["theta"] = df["theta_main"] + df["theta_aux"] / 60
df["theta_ref"] = df["theta_main_ref"] + df["theta_aux_ref"] / 60
df["theta_obs"] = np.abs(df["theta"] - df["theta_ref"])

# =========================
# Experiment A
# =========================
LAMBDA_REF = {3: lambda_ref_3, 4: lambda_ref_4, 5: lambda_ref_5, 6: lambda_ref_6}
CAPTION_BY_M = {1: "第一級繞射", 2: "第二級繞射", 3: "第三級繞射"}

results_A = {}
R_list = []

for m_val in sorted(df[df["subexp_id"] == "A"]["m"].unique()):
    df_m = df[(df["subexp_id"] == "A") & (df["m"] == m_val)]
    grouped = df_m.groupby("ni")

    table_rows = []

    for ni, g in grouped:
        theta_list = g["theta_obs"].values
        theta_str = ", ".join([f"{v:.2f}" for v in theta_list])
        theta_rad = np.deg2rad(g["theta_obs"])
        lambda_val = np.mean(np.sin(theta_rad)) * d / m_val

        lambda_ref_val = LAMBDA_REF.get(ni, np.nan)
        rel_err = (lambda_val - lambda_ref_val) / lambda_ref_val * 100
        R_val = 1 / (lambda_val * (1/4 - 1/(ni**2)))
        R_list.append(R_val)

        table_rows.append({
            "ni": ni,
            "thetaobs": theta_str,
            "lambda": lambda_val,
            "lambdaref": lambda_ref_val,
            "relerrpercent": rel_err,
            "R": R_val,
        })

    results_A[m_val] = pd.DataFrame(table_rows)

# R average and uncertainty
R_array = np.array(R_list)
R_avg = np.mean(R_array) * 1e9
u_R = np.std(R_array, ddof=1) / np.sqrt(len(R_array)) * 1e9
R_rel_err = (R_avg - R_ref) / R_ref * 100

# Planck constant from R
h_est = np.cbrt((e**4 * mu) / (8 * eps0**2 * c * R_avg))
u_h = h_est * (1/3) * (u_R / R_avg)
h_rel_err = (h_est - h_ref) / h_ref * 100

# A tables per m
for m_val, table in results_A.items():
    table_out = pd.DataFrame({
        "$n_{i}$": table["ni"].map(lambda x: f"${x:.4g}$"),
        "繞射角$\\theta_{\\text{obs}}$ (\\unit{\\degree})": table["thetaobs"],
        "$\\lambda_{\\text{obs}}$ (\\unit{\\nm})": table["lambda"].map(lambda x: f"${x:.4g}$"),
        "$\\lambda_{\\text{ref}}$ (\\unit{\\nm})": table["lambdaref"].map(lambda x: f"${x:.4g}$"),
        "Error": table["relerrpercent"].map(percent),
        "Rydberg constant (\\unit{\\per\\nm})": table["R"].map(lambda x: f"${x:.6g}$"),
    })
    df_to_latex(table_out, path=REPORT, caption=CAPTION_BY_M.get(m_val, f"第{m_val}級繞射"))

# =========================
# linear regression (m=1)
# =========================
df_m1 = results_A[1]
ni_vals = df_m1["ni"].values
lambda_vals = df_m1["lambda"].values

x1 = 1 / (ni_vals**2)
y1 = 1 / lambda_vals

popt, pcov = curve_fit(linear_model, x1, y1)
s, b = popt
u_s, u_b = np.sqrt(np.diag(pcov))

ss_tot = np.sum((y1 - np.mean(y1))**2)
R2 = 1 - np.sum((y1 - linear_model(x1, *popt))**2) / ss_tot

a = -b / s
u_a = np.sqrt((u_b/s)**2 + (b*u_s/s**2)**2 + 2*b/s**3*pcov[0, 1])

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(x1, y1, label="Data Point")
x_fit = np.linspace(min(x1), max(x1), 100)
if b >= 0:
    fit_label = f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {R2:.4f})'
else:
    fit_label = f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {R2:.4f})'
plt.plot(x_fit, linear_model(x_fit, *popt), color='black', linestyle='--', label=fit_label)
plt.legend()
plt.xlabel(r"$1/n_i^2$")
plt.ylabel(r"$1/\lambda$ (nm$^{-1}$)")
plt.savefig(FIGURES_DIR / "inv_ni2_1.png", dpi=300)
plt.close()

neg_1e9_s = -1e9 * s
u_neg_1e9_s = 1e9 * u_s
err_s = (neg_1e9_s - R_ref) / R_ref * 100

four_1e9_b = 4e9 * b
u_four_1e9_b = 4e9 * u_b
err_b = (four_1e9_b - R_ref) / R_ref * 100

err_a = (a - 0.25) / 0.25 * 100

x2 = 1/4 - 1/(ni_vals**2)
popt2, pcov2 = curve_fit(linear_model, x2, y1)
s2, b2 = popt2
R2_2 = 1 - np.sum((y1 - linear_model(x2, *popt2))**2) / ss_tot

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(x2, y1, label="Data Point")
x_fit2 = np.linspace(min(x2), max(x2), 100)
if b2 >= 0:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x + {b2:.4g} (R² = {R2_2:.4f})'
else:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x - {abs(b2):.4g} (R² = {R2_2:.4f})'
plt.plot(x_fit2, linear_model(x_fit2, *popt2), color='black', linestyle='--', label=fit_label2)
plt.legend()
plt.xlabel(r"$1/4 - 1/n_i^2$")
plt.ylabel(r"$1/\lambda$ (nm$^{-1}$)")
plt.savefig(FIGURES_DIR / "transition_1.png", dpi=300)
plt.close()

fit_compare_m1 = pd.DataFrame({
    "Quantity": [
        "$-s$ (\\unit{\\per\\m})",
        "$4b$ (\\unit{\\per\\m})",
        "x-intercept",
    ],
    "Value": [
        pm(neg_1e9_s, u_neg_1e9_s),
        pm(four_1e9_b, u_four_1e9_b),
        pm(a, u_a),
    ],
    "Error": [
        percent(err_s),
        percent(err_b),
        percent(err_a),
    ],
})
df_to_latex(fit_compare_m1, path=REPORT, caption="fitting result for $m=1$")

# =========================
# linear regression (m=2)
# =========================
df_m2 = results_A[2]
ni_vals = df_m2["ni"].values
lambda_vals = df_m2["lambda"].values

x1 = 1 / (ni_vals**2)
y1 = 1 / lambda_vals

popt, pcov = curve_fit(linear_model, x1, y1)
s, b = popt
u_s, u_b = np.sqrt(np.diag(pcov))

ss_tot = np.sum((y1 - np.mean(y1))**2)
R2 = 1 - np.sum((y1 - linear_model(x1, *popt))**2) / ss_tot

a = -b / s
u_a = np.sqrt((u_b/s)**2 + (b*u_s/s**2)**2 + 2*b/s**3*pcov[0, 1])

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(x1, y1, label="Data Point")
x_fit = np.linspace(min(x1), max(x1), 100)
if b >= 0:
    fit_label = f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {R2:.4f})'
else:
    fit_label = f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {R2:.4f})'
plt.plot(x_fit, linear_model(x_fit, *popt), color='black', linestyle='--', label=fit_label)
plt.legend()
plt.xlabel(r"$1/n_i^2$")
plt.ylabel(r"$1/\lambda$ (nm$^{-1}$)")
plt.savefig(FIGURES_DIR / "inv_ni2_2.png", dpi=300)
plt.close()

neg_1e9_s = -1e9 * s
u_neg_1e9_s = 1e9 * u_s
err_s = (neg_1e9_s - R_ref) / R_ref * 100

four_1e9_b = 4e9 * b
u_four_1e9_b = 4e9 * u_b
err_b = (four_1e9_b - R_ref) / R_ref * 100

err_a = (a - 0.25) / 0.25 * 100

x2 = 1/4 - 1/(ni_vals**2)
popt2, pcov2 = curve_fit(linear_model, x2, y1)
s2, b2 = popt2
R2_2 = 1 - np.sum((y1 - linear_model(x2, *popt2))**2) / ss_tot

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(x2, y1, label="Data Point")
x_fit2 = np.linspace(min(x2), max(x2), 100)
if b2 >= 0:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x + {b2:.4g} (R² = {R2_2:.4f})'
else:
    fit_label2 = f'Best-fit line: y = {s2:.4g} x - {abs(b2):.4g} (R² = {R2_2:.4f})'
plt.plot(x_fit2, linear_model(x_fit2, *popt2), color='black', linestyle='--', label=fit_label2)
plt.legend()
plt.xlabel(r"$1/4 - 1/n_i^2$")
plt.ylabel(r"$1/\lambda$ (nm$^{-1}$)")
plt.savefig(FIGURES_DIR / "transition_2.png", dpi=300)
plt.close()

fit_compare_m2 = pd.DataFrame({
    "Quantity": [
        "$-s$ (\\unit{\\per\\m})",
        "$4b$ (\\unit{\\per\\m})",
        "x-intercept",
    ],
    "Value": [
        pm(neg_1e9_s, u_neg_1e9_s),
        pm(four_1e9_b, u_four_1e9_b),
        pm(a, u_a),
    ],
    "Error": [
        percent(err_s),
        percent(err_b),
        percent(err_a),
    ],
})
df_to_latex(fit_compare_m2, path=REPORT, caption="fitting result for $m=2$")

# =========================
# Experiment B
# =========================
df_B = df[df["subexp_id"] == "B"].copy()
df_B["n"] = np.sin((alpha + np.deg2rad(df_B["theta_obs"])) / 2) / np.sin(alpha/2)

xB = results_A[1]["lambdaref"].values
yB = df_B["n"].values

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(xB, yB, label="Data Point")
plt.legend()
plt.xlabel(r"$\lambda$ (nm)")
plt.ylabel("$n$")
plt.savefig(FIGURES_DIR / "n_vs_lambda.png", dpi=300)
plt.close()

summary_constants = pd.DataFrame({
    "Quantity": [
        "$R_{\\text{est}}$ (\\unit{\\per\\m})",
        "$R_{\\text{theo}}=\\dfrac{\\mu e^{4}}{8\\varepsilon_{0}^{2}h_{\\text{ref}}^{3}c}$ (\\unit{\\per\\m})",
        "Error",
        "$h_{\\text{est}}$ (\\unit{\\J\\s})",
        "$h_{\\text{ref}}$ (\\unit{\\J\\s})",
        "Error",
    ],
    "Value": [
        pm(R_avg, u_R),
        f"${R_ref:.4g}$",
        percent(R_rel_err),
        pm(h_est, u_h),
        f"${h_ref:.4g}$",
        percent(h_rel_err),
    ],
})
df_to_latex(summary_constants, path=REPORT, caption="物理常數的計算結果")

df_B_out = pd.DataFrame({
    "$n_{i}$": df_B["ni"].map(lambda x: f"${x:.4g}$"),
    "最小偏折角$D$": df_B["theta_obs"].map(lambda x: f"${x:.4g}$"),
    "折射率$n$": df_B["n"].map(lambda x: f"${x:.4g}$"),
})
df_to_latex(df_B_out, path=REPORT, caption="不同譜線對應的折射率")
