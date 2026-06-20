import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, linear_model

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 14" / "figures"
REPORT = BASE_DIR / "lab 14" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

df = pd.read_excel(BASE_DIR / "lab 14" / "data.xlsx")

g = 9.80665
I_const_A = 3.04
B_ref = 0.09677 - 0.00064

D = 0.0177
a_param = 0.2115
b_param = 1.0065
R = 0.0031
L = 0.279
mu0 = 4 * np.pi * 10**(-7)

I_ref_B = I_const_A

df["mass"] = df["mass"] / 1000

# =========================
# Sub-experiment A: F vs length
# =========================
df_A = df[df["subexp"] == "A"].copy()
df_A["F"] = df_A["mass"] * g
df_A["outer_length"] = df_A["outer_length"] / 1000
df_A["inner_length"] = df_A["inner_length"] / 1000
df_A["average_length"] = (df_A["outer_length"] + df_A["inner_length"]) / 2

popt_outer, pcov_outer = curve_fit(linear_model, df_A["outer_length"], df_A["F"])
popt_inner, pcov_inner = curve_fit(linear_model, df_A["inner_length"], df_A["F"])
popt_avg,   pcov_avg   = curve_fit(linear_model, df_A["average_length"], df_A["F"])

s_outer, b_outer = popt_outer
s_inner, b_inner = popt_inner
s_avg,   b_avg   = popt_avg

u_s_outer = np.sqrt(pcov_outer[0, 0])
u_s_inner = np.sqrt(pcov_inner[0, 0])
u_s_avg   = np.sqrt(pcov_avg[0, 0])


def plot_fit(x, y, s, b, xlabel, ylabel, figname):
    x = np.asarray(x)
    y = np.asarray(y)
    y_predicted = s * x + b
    residuals = y - y_predicted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = linear_model(x_line, s, b)
    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.scatter(x, y, label="Data Point")
    if b < 0:
        plt.plot(x_line, y_line, color='black', linestyle='--',
                 label=f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {r_squared:.4f})')
    else:
        plt.plot(x_line, y_line, color='black', linestyle='--',
                 label=f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {r_squared:.4f})')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(FIGURES_DIR / (figname + ".png"), dpi=300)
    plt.close()


plot_fit(df_A["outer_length"],   df_A["F"], s_outer, b_outer,
         "outer length (m)",   "force (N)", "outer length")
plot_fit(df_A["inner_length"],   df_A["F"], s_inner, b_inner,
         "inner length (m)",   "force (N)", "inner length")
plot_fit(df_A["average_length"], df_A["F"], s_avg,   b_avg,
         "average length (m)", "force (N)", "average length")

out_x = np.asarray(df_A["outer_length"])
inn_x = np.asarray(df_A["inner_length"])
avg_x = np.asarray(df_A["average_length"])
y_F   = np.asarray(df_A["F"])

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.scatter(out_x, y_F, label="Using outer")
plt.scatter(inn_x, y_F, label="Using inner")
plt.scatter(avg_x, y_F, label="Using average")
plt.xlabel("length (m)")
plt.ylabel("force (N)")
plt.legend()
plt.savefig(FIGURES_DIR / "lengths.png", dpi=300)
plt.close()

B_outer = s_outer / I_const_A
B_inner = s_inner / I_const_A
B_avg   = s_avg   / I_const_A

u_B_outer = u_s_outer / I_const_A
u_B_inner = u_s_inner / I_const_A
u_B_avg   = u_s_avg   / I_const_A

err_outer = (B_outer - B_ref) / B_ref * 100
err_inner = (B_inner - B_ref) / B_ref * 100
err_avg   = (B_avg   - B_ref) / B_ref * 100

# =========================
# Sub-experiment B: F vs B (varying magnet count)
# =========================
df_B = df[df["subexp"] == "B"].copy()
df_B["F"] = df_B["mass"] * g

popt_B, pcov_B = curve_fit(linear_model, df_B["B_values"], df_B["F"])
s_mag, b_mag = popt_B
u_s_mag = np.sqrt(pcov_B[0, 0])
u_b_mag = np.sqrt(pcov_B[1, 1])

current   = s_mag / 0.038
u_current = u_s_mag / 0.038
err_i     = (current - I_ref_B) / I_ref_B * 100

plot_fit(df_B["B_values"], df_B["F"], s_mag, b_mag,
         "B (T)", "force (N)", "magnetic field")

# =========================
# Sub-experiment C: F vs I² (current balance, measure µ)
# =========================
df_C = df[df["subexp"] == "C"].copy()
df_C["F"]  = df_C["mass"] * g
df_C["I2"] = df_C["I_power_supply"] ** 2

popt_C, pcov_C = curve_fit(linear_model, df_C["I2"], df_C["F"])
s_C, b_C   = popt_C
var_s_C    = pcov_C[0, 0]
u_s_C      = np.sqrt(var_s_C)
var_b_C    = pcov_C[1, 1]
cov_sb_C   = pcov_C[0, 1]

plot_fit(df_C["I2"], df_C["F"], s_C, b_C,
         "I$^2$ ($\\mathrm{A}^2$)", "force (N)", "Squared current")

d               = (a_param * D / (2 * b_param)) + R
K_mu            = 2 * np.pi * d / L
mu_est          = K_mu * s_C
u_mu_est        = K_mu * u_s_C
mu_error_percent = (mu_est - mu0) / mu0 * 100

# =========================
# Sub-experiment D: predict unknown mass
# =========================
df_D = df[df["subexp"] == "D"].copy()
df_D["I2"] = df_D["I_power_supply"] ** 2
x = df_D["I2"]

F_pred     = linear_model(x, s_C, b_C)
mass_prime = F_pred / g

var_F        = x**2 * var_s_C + var_b_C + 2 * x * cov_sb_C
u_F          = np.sqrt(var_F)
u_mass_prime = u_F / g

df_D["F_pred"]            = F_pred
df_D["u_F"]               = u_F
df_D["mass_prime"]        = mass_prime
df_D["u_mass_prime"]      = u_mass_prime
df_D["percent_error_mass"] = (df_D["mass_prime"] - df_D["mass"]) / df_D["mass"] * 100

# =========================
# Output tables
# =========================

df_to_latex(
    pd.DataFrame({
        "迴路板編號": df_A["board_number"].map(str),
        "outer length (\\unit{\\mm})":   df_A["outer_length"].map(lambda x: f"${1000*x:.1f}$"),
        "inner length (\\unit{\\mm})":   df_A["inner_length"].map(lambda x: f"${1000*x:.1f}$"),
        "average length (\\unit{\\mm})": df_A["average_length"].map(lambda x: f"${1000*x:.1f}$"),
        "force (gw)":          df_A["mass"].map(lambda x: f"${1000*x:.4g}$"),
        "force (\\unit{\\N})": df_A["F"].map(lambda x: f"${x:.3g}$"),
    }),
    caption="Measured lengths and Lorentz force with $I_{\\text{obs}}=\\qty{3.04}{\\A}$",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "slope $s$",
            "intercept",
            "$B=s/I_{\\text{obs}}$ (\\unit{\\T})",
            "error",
        ],
        "Outer": [
            pm(s_outer, u_s_outer),
            pm(b_outer, np.sqrt(pcov_outer[1, 1])),
            pm(B_outer, u_B_outer),
            percent(err_outer),
        ],
        "Average": [
            pm(s_avg, u_s_avg),
            pm(b_avg, np.sqrt(pcov_avg[1, 1])),
            pm(B_avg, u_B_avg),
            percent(err_avg),
        ],
        "Inner": [
            pm(s_inner, u_s_inner),
            pm(b_inner, np.sqrt(pcov_inner[1, 1])),
            pm(B_inner, u_B_inner),
            percent(err_inner),
        ],
    }),
    caption="Linear fit results with $B_{\\text{obs}} = \\qty{96.13}{\\mT}$",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "magnet number":       df_B["magnet_number"].map(lambda x: f"${x:.0f}$"),
        "B (\\unit{\\mT})":  df_B["B_values"].map(lambda x: f"${1000*x:.4g}$"),
        "force (gw)":          df_B["mass"].map(lambda x: f"${1000*x:.2f}$"),
        "force (\\unit{\\N})": df_B["F"].map(lambda x: f"${x:.4g}$"),
    }),
    caption="改變磁鐵數量，測量磁場與洛倫茲力的大小。使用38號板",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "slope $s$",
            "intercept",
            "$I=s/L$ (\\unit{\\A})",
            "$I_{\\text{obs}}$ (\\unit{\\A})",
            "error",
        ],
        "Value": [
            pm(s_mag, u_s_mag),
            pm(b_mag, u_b_mag),
            pm(current, u_current),
            "$3.04$",
            percent(err_i),
        ],
    }),
    caption="線性擬合的結果。$L= \\qty{38.0}{\\mm}$ 是38號板的內徑，",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "mass (\\unit{\\g})":       df_C["mass"].map(lambda x: f"${1000*x:.4f}$"),
        "force (\\qty{e-4}{\\N})":  df_C["F"].map(lambda x: f"${1e4*x:.4g}$"),
        "current (\\unit{\\A})":    df_C["I_power_supply"].map(lambda x: f"${x:.3f}$"),
    }),
    caption="Measured mass and calculated force versus power-supply current",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "slope $s$",
            "intercept",
            "反射鏡到動導線的距離 $a$ (\\unit{\\cm})",
            "反射鏡到標尺的距離 $b$ (\\unit{\\cm})",
            "光點位移 $D$ (\\unit{\\cm})",
            "動導線直徑 $2r_{0}$ (\\unit{\\mm})",
            "兩導線中心的距離 $d=2r_{0}+aD/2b$ (\\unit{\\m})",
            "動導線長度 $L$ (\\unit{\\cm})",
            "$\\mu_{\\text{est}}=2\\pi sd/L$ (\\unit{\\N/\\square\\A})",
            "$\\mu_0$ (\\unit{\\N/\\square\\A})",
            "error",
        ],
        "Value": [
            pm(s_C, u_s_C),
            pm(b_C, np.sqrt(var_b_C)),
            f"${a_param * 100:.2f}$",
            f"${b_param * 100:.2f}$",
            f"${D * 100:.2f}$",
            f"${R * 1000:.2f}$",
            f"${d:.3g}$",
            f"${L * 100:.2f}$",
            pm(mu_est, u_mu_est),
            "$1.257\\times 10^{-6}$",
            percent(mu_error_percent),
        ],
    }),
    caption="Derived parameters from the linear fit and measurement",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "current (\\unit{\\A})":          df_D["I_power_supply"].map(lambda x: f"${x:.4g}$"),
        "force (\\unit{\\N})":            [pm(v, u) for v, u in zip(df_D["F_pred"], df_D["u_F"])],
        "predicted mass (\\unit{\\kg})":  [
            pm(v, u) for v, u in zip(df_D["mass_prime"], df_D["u_mass_prime"])
        ],
        "measured mass (\\unit{\\kg})":   df_D["mass"].map(lambda x: f"${x:.4g}$"),
        "error":                          df_D["percent_error_mass"].map(percent),
    }),
    caption="Predicted mass from fitted model and comparison with measured mass",
    path=REPORT,
)
