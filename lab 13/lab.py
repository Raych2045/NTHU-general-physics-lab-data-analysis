import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
from scipy.odr import ODR, Model, RealData

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, linear_model

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 13" / "figures"
REPORT = BASE_DIR / "lab 13" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

df = pd.read_excel(BASE_DIR / "lab 13" / "data.xlsx")
df["I_Galvanometer"] = df["I_Galvanometer"] * 1e-6
df["I_multimeter"] = df["I_multimeter"] * 1e-3

u_IG = 2/np.sqrt(12) * 1e-6

# =========================
# Sub-experiment A: calibrate galvanometer
# =========================
df_A = df[df["subexp_id"] == "A"].copy()

df_A["I_T"] = df_A["V_multimeter"] / df_A["R_load"]

u_R_ext = df_A["R_load"] * df_A["uncertainty_R_load(%)"] / 100
df_A["u_I_T"] = df_A["I_T"] * (u_R_ext / df_A["R_load"])

x = df_A["I_Galvanometer"].to_numpy()
y = df_A["I_T"].to_numpy()
yerr = df_A["u_I_T"].to_numpy()

def linear_model_forA(B, x):
    return B[0]*x + B[1]

model = Model(linear_model_forA)
data = RealData(x, y, sx=u_IG, sy=yerr)

odr = ODR(data, model, beta0=[1, 0])
output = odr.run()

a, b = output.beta
pcov_A = output.cov_beta * output.res_var

var_a, var_b = np.diag(pcov_A)
cov_ab = pcov_A[0, 1]

y_calculated = linear_model_forA(output.beta, x)
residuals = y - y_calculated
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

x_fit = np.linspace(x.min(), x.max(), 200)
y_fit = linear_model_forA(output.beta, x_fit)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(x, y, xerr=u_IG, yerr=yerr, fmt="o", capsize=3, label="Data Point")
if b < 0:
    plt.plot(x_fit, y_fit, color='black', linestyle="--",
             label=f'Best-fit line: y = {a:.4g} x - {abs(b):.4g} (R² = {r_squared:.4f})')
else:
    plt.plot(x_fit, y_fit, color='black', linestyle="--",
             label=f'Best-fit line: y = {a:.4g} x + {b:.4g} (R² = {r_squared:.4f})')
plt.xlabel("I_Galvanometer (A)")
plt.ylabel("I_T (A)")
plt.legend()
plt.savefig(FIGURES_DIR / "I_T-I_G.png", dpi=300)
plt.close()

# =========================
# Sub-experiment B: measure galvanometer resistance
# =========================
df_B = df[df["subexp_id"] == "B"].copy()
I_G = df_B["I_Galvanometer"]
df_B["I_T"] = linear_model(I_G, a, b)
df_B["u_I_T"] = np.sqrt(
    (I_G**2) * var_a +
    var_b +
    2 * I_G * cov_ab +
    a**2 * u_IG**2
)

x = df_B["V_multimeter"].to_numpy()
y = df_B["I_T"].to_numpy()
yerr = df_B["u_I_T"].to_numpy()

popt_B, pcov_B = curve_fit(
    linear_model,
    x, y,
    sigma=yerr,
    absolute_sigma=True
)

RG_reciprocal, intercept_B = popt_B
u_RG_reciprocal = np.sqrt(pcov_B[0, 0])
R_G = 1/RG_reciprocal
u_R_G = u_RG_reciprocal / RG_reciprocal**2

y_calculated = RG_reciprocal * x + intercept_B
residuals = y - y_calculated
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

x_fit = np.linspace(x.min(), x.max(), 200)
y_fit = linear_model(x_fit, RG_reciprocal, intercept_B)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, label="Data Point")
if intercept_B < 0:
    plt.plot(x_fit, y_fit, color='black', linestyle="--",
             label=f'Best-fit line: y = {RG_reciprocal:.4g} x - {abs(intercept_B):.4g} (R² = {r_squared:.4f})')
else:
    plt.plot(x_fit, y_fit, color='black', linestyle="--",
             label=f'Best-fit line: y = {RG_reciprocal:.4g} x + {intercept_B:.4g} (R² = {r_squared:.4f})')
plt.xlabel("V_Multimeter (V)")
plt.ylabel("I_T (A)")
plt.legend()
plt.savefig(FIGURES_DIR / "V_PS-I_T.png", dpi=300)
plt.close()

# =========================
# Sub-experiment C: self-made ammeter
# =========================
df_C = df[df["subexp_id"] == "C"].copy()

I_G = df_C["I_Galvanometer"].to_numpy()

df_C["I_T"] = linear_model(I_G, a, b)
df_C["u_I_T"] = np.sqrt(
    (I_G**2) * var_a +
    var_b +
    2 * I_G * cov_ab +
    a**2 * u_IG**2
)

R_used = df_C["R_used"]
u_R_used = [1.65, np.sqrt(15**2+1+1), np.sqrt(15**2+1)] * df_C["uncertainty_R_used(%)"] / 100

df_C["I"] = df_C["I_T"] * (1 + R_G / R_used)
df_C["u_I"] = np.sqrt(
    ((1 + R_G / R_used) * df_C["u_I_T"])**2 +
    ((df_C["I_T"] / R_used) * u_R_G)**2 +
    ((df_C["I_T"] * R_G / R_used**2) * u_R_used)**2
)
df_C["percent_error_I"] = (
    (df_C["I"] - df_C["I_multimeter"]) / df_C["I_multimeter"] * 100
)

# =========================
# Sub-experiment D: self-made voltmeter
# =========================
df_D = df[df["subexp_id"] == "D"].copy()

I_G = df_D["I_Galvanometer"].to_numpy()

df_D["I_T"] = linear_model(I_G, a, b)
df_D["u_I_T"] = np.sqrt(
    (I_G**2) * var_a +
    var_b +
    2 * I_G * cov_ab +
    a**2 * u_IG**2
)

R_used = df_D["R_used"]
u_R_used = R_used * df_D["uncertainty_R_used(%)"] / 100

df_D["V"] = df_D["I_T"] * (R_G + R_used)
df_D["u_V"] = np.sqrt(
    ((R_G + R_used) * df_D["u_I_T"])**2 +
    (df_D["I_T"] * u_R_G)**2 +
    (df_D["I_T"] * u_R_used)**2
)
df_D["percent_error_V"] = (
    (df_D["V"] - df_D["V_multimeter"]) / df_D["V_multimeter"] * 100
)

# =========================
# Sub-experiment E: self-made ohmmeter
# =========================
df_E = df[df["subexp_id"] == "E"].copy()

I_G = df_E["I_Galvanometer"].to_numpy()

df_E["I_T"] = linear_model(I_G, a, b)
df_E["u_I_T"] = np.sqrt(
    (I_G**2) * var_a +
    var_b +
    2 * I_G * cov_ab +
    a**2 * u_IG**2
)

R_used = df_E["R_used"]
u_R_used = R_used * df_E["uncertainty_R_used(%)"] / 100.0

df_E["R_prime"] = (
    (0.00005 - df_E["I_T"]) * (R_G + R_used) / df_E["I_T"]
)

dR_dRG = (0.00005 - df_E["I_T"]) / df_E["I_T"]
dR_dIG = -0.00005 * (R_G + R_used) / (df_E["I_T"]**2)

df_E["u_R_prime"] = np.sqrt(
    (dR_dRG * u_R_G)**2 +
    (dR_dRG * u_R_used)**2 +
    (dR_dIG * df_E["u_I_T"])**2
)

df_E["percent_error_R"] = (
    (df_E["R_prime"] - df_E["R_load_multimeter"]) / df_E["R_load_multimeter"] * 100
)

# =========================
# Output tables
# =========================

df_to_latex(
    pd.DataFrame({
        "$V_{\\text{multi}}$ (\\unit{\\V})":
            df_A["V_multimeter"].map(lambda x: f"${x:.2f}$"),
        "$I_G$ (\\unit{\\uA})":
            df_A["I_Galvanometer"].map(lambda x: f"${1e6*x:.1f}$"),
        "$I_T$ (\\unit{\\uA})": [
            pm(v, u) for v, u in zip(1e6*df_A["I_T"], 1e6*df_A["u_I_T"])
        ],
    }),
    caption=(
        "Measured and calculated current, using "
        "$R_{\\text{load}}=\\qty{390}{\\kohm}$ with 5\\% tolerance"
    ),
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": ["$a$", "$b$"],
        "Value": [
            pm(a, np.sqrt(var_a)),
            pm(b, np.sqrt(var_b)),
        ],
    }),
    caption="Linear fit parameters for $I_T = a I_G + b$",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "$V_{\\text{multi}}$ (\\unit{\\mV})":
            df_B["V_multimeter"].map(lambda x: f"${x*1e3:.1f}$"),
        "$I_G$ (\\unit{\\uA})":
            df_B["I_Galvanometer"].map(lambda x: f"${1e6*x:.1f}$"),
        "$I_T$ (\\unit{\\uA})": [
            pm(v, u) for v, u in zip(1e6*df_B["I_T"], 1e6*df_B["u_I_T"])
        ],
    }),
    caption="Measured voltage and corrected current",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "$a$",
            "$b$",
            "$R_G=1/a$ (\\unit{\\kohm})",
        ],
        "Value": [
            pm(RG_reciprocal, u_RG_reciprocal),
            pm(intercept_B, np.sqrt(pcov_B[1, 1])),
            pm(R_G/1e3, u_R_G/1e3),
        ],
    }),
    caption="Calculating galvanometer resistance $R_G$",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "$R_{\\text{used}}$ (\\unit{\\ohm})":
            df_C["R_used"].map(lambda x: f"${x:.4g}$"),
        "$R_{\\text{load}}$ (\\unit{\\ohm})":
            df_C["R_load"].map(lambda x: f"${x:.4g}$"),
        "$V_{PS}$ (\\unit{\\V})":
            df_C["V_power_supply"].map(lambda x: f"${x:.2f}$"),
        "$I_G$ (\\unit{\\uA})":
            df_C["I_Galvanometer"].map(lambda x: f"${1e6*x:.4g}$"),
        "$I_T$ (\\unit{\\uA})": [
            pm(v, u) for v, u in zip(1e6*df_C["I_T"], 1e6*df_C["u_I_T"])
        ],
        "$I_{\\text{tot}}=I_{T}(1+R_{G}/R_{\\text{used}})$ (\\unit{\\A})": [
            pm(v, u) for v, u in zip(df_C["I"], df_C["u_I"])
        ],
        "$I_{\\text{multi}}$ (\\unit{\\A})":
            df_C["I_multimeter"].map(lambda x: f"${x:.4g}$"),
        "Error":
            df_C["percent_error_I"].map(percent),
    }),
    caption="Current comparison between self-made ammeter and multimeter",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "$R_{\\text{used}}$ (\\unit{\\kohm})":
            df_D["R_used"].map(lambda x: f"${x/1e3:.4g}$"),
        "$R_{\\text{load}}$ (\\unit{\\ohm})":
            df_D["R_load"].map(lambda x: f"${x:.4g}$"),
        "$V_{PS}$ (\\unit{\\V})":
            df_D["V_power_supply"].map(lambda x: f"${x:.2f}$"),
        "$I_G$ (\\unit{\\uA})":
            df_D["I_Galvanometer"].map(lambda x: f"${1e6*x:.4g}$"),
        "$I_T$ (\\unit{\\uA})": [
            pm(v, u) for v, u in zip(1e6*df_D["I_T"], 1e6*df_D["u_I_T"])
        ],
        "$V_{\\text{ab}}=I_{T}(R_{G}+R_{\\text{used}})$ (\\unit{\\V})": [
            pm(v, u) for v, u in zip(df_D["V"], df_D["u_V"])
        ],
        "$V_{\\text{multi}}$ (\\unit{\\V})":
            df_D["V_multimeter"].map(lambda x: f"${x:.4g}$"),
        "Error":
            df_D["percent_error_V"].map(percent),
    }),
    caption="Voltage comparison between self-made voltmeter and multimeter",
    path=REPORT,
)

for Vps, g in df_E.groupby("V_multimeter"):
    r_used_kohm = g["R_used"].iloc[0] / 1e3
    df_to_latex(
        pd.DataFrame({
            "$R_{\\text{load}}$ (\\unit{\\kohm})":
                g["R_load"].map(lambda x: f"${x/1e3:.4g}$"),
            "$I_G$ (\\unit{\\uA})":
                g["I_Galvanometer"].map(lambda x: f"${1e6*x:.4g}$"),
            "$I_T$ (\\unit{\\uA})": [
                pm(v, u) for v, u in zip(1e6*g["I_T"], 1e6*g["u_I_T"])
            ],
            "$R_{\\text{load,est}}$ (\\unit{\\kohm})": [
                pm(v, u) for v, u in zip(g["R_prime"]/1e3, g["u_R_prime"]/1e3)
            ],
            "$R_{\\text{load,multi}}$ (\\unit{\\kohm})":
                g["R_load_multimeter"].map(lambda x: f"${x/1e3:.4g}$"),
            "Error":
                g["percent_error_R"].map(percent),
        }),
        caption=(
            f"Resistance comparison between self-made ohmmeter and multimeter "
            f"with $V_{{PS}}={Vps:.4g}\\,\\mathrm{{V}}$ and "
            f"$R_{{\\text{{used}}}}=\\qty{{{r_used_kohm:.4g}}}{{\\kohm}}$"
        ),
        path=REPORT,
    )
