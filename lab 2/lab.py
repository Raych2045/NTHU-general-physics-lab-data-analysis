import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty, linear_model, compute_r_squared

G_CGS = 980.665               # cm/s^2, standard gravity (this lab works in cm-g-s units)
M_CART_REF_G = 382.70       # g, cart mass measured directly by an electronic balance ("other method")
SLOPE_RULER_LENGTH_CM = 113.92   # cm, ruler-measured incline length
SLOPE_HEIGHT_DIFF_CM = 4.55      # cm, ruler-measured height difference between the incline's two ends
DYNE_PER_N = 1e5             # 1 N = 1e5 dyne (g*cm/s^2), since masses are in g and accelerations in cm/s^2
N_FIT_A_PLUS_LINEAR = 15     # only the first 15 (smallest-m) a_+ points stay linear
ACCEL_RES_CM_S2 = 0.01       # cm/s^2, photogate timer's acceleration-mode reading resolution

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 2" / "figures"
REPORT = BASE_DIR / "lab 2" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

DATA_PATH = BASE_DIR / "lab 2" / "data.xlsx"

df_incline = pd.read_excel(DATA_PATH, sheet_name="incline")
df_horizontal = pd.read_excel(DATA_PATH, sheet_name="horizontal")

sin_theta_other = SLOPE_HEIGHT_DIFF_CM / SLOPE_RULER_LENGTH_CM
theta_other_deg = np.degrees(np.arcsin(sin_theta_other))


def plot_scatter_fit(x, y, y_err, x_fit_curve, y_fit_curve, fit_label, xlabel, ylabel, title, filename):
    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    if y_err is not None:
        plt.errorbar(x, y, yerr=y_err, fmt="o", capsize=3, label="Data Point")
    else:
        plt.scatter(x, y, label="Data Point")
    if x_fit_curve is not None:
        plt.plot(x_fit_curve, y_fit_curve, color='black', linestyle='--', label=fit_label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(FIGURES_DIR / filename, dpi=300)
    plt.close()


def accel_stats(arr):
    mean, uA = typeA_uncertainty(arr)
    uB = ACCEL_RES_CM_S2 / np.sqrt(12)
    uC = np.sqrt(uA**2 + uB**2)
    return mean, uC


def linear_fit_with_xintercept(x, y, sigma):
    popt, pcov = curve_fit(linear_model, x, y, sigma=sigma, absolute_sigma=True)
    s, b = popt
    var_s, var_b = pcov[0, 0], pcov[1, 1]
    cov_sb = pcov[0, 1]
    u_s, u_b = np.sqrt(var_s), np.sqrt(var_b)
    X = -b / s
    u_X = np.sqrt((b / s**2)**2 * var_s + (1 / s)**2 * var_b - 2 * (b / s**2) * (1 / s) * cov_sb)
    r2 = compute_r_squared(x, s, b, y)
    return s, u_s, b, u_b, X, u_X, r2


# =========================
# Tables 1 & 2: incline-motion acceleration vs. hanging mass
# =========================
trials = df_incline["acceleration_cm_s2"].apply(lambda s: np.array([float(v) for v in str(s).split(",")]))
stats = trials.apply(accel_stats)
mean_a = stats.apply(lambda t: t[0]).to_numpy()
u_a = stats.apply(lambda t: t[1]).to_numpy()
mass = df_incline["mass_g"].to_numpy()
is_minus = (df_incline["group"] == "a_minus").to_numpy()
is_plus = (df_incline["group"] == "a_plus").to_numpy()

df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": [f"${v:.2f}$" for v in mass[is_minus]],
        "Acceleration $a$ (\\unit{\\cm\\per\\s\\squared})": [pm(m_, u_) for m_, u_ in zip(mean_a[is_minus], u_a[is_minus])],
    }),
    caption="Incline motion -- acceleration vs. hanging mass ($a_{-}$ group, including the zero-crossing points)",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": [f"${v:.2f}$" for v in mass[is_plus]],
        "Acceleration $a$ (\\unit{\\cm\\per\\s\\squared})": [pm(m_, u_) for m_, u_ in zip(mean_a[is_plus], u_a[is_plus])],
    }),
    caption="Incline motion -- acceleration vs. hanging mass ($a_{+}$ group)",
    path=REPORT,
)

# Overview: a_- values negated, a_+ values as-is, combined into one a-m graph with error bars
m_overview = np.concatenate([mass[is_minus], mass[is_plus]])
a_overview = np.concatenate([-mean_a[is_minus], mean_a[is_plus]])
u_overview = np.concatenate([u_a[is_minus], u_a[is_plus]])

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(m_overview, a_overview, yerr=u_overview, fmt="o", capsize=3, label="Data Point")
plt.xlabel("Hanging mass $m$ (g)")
plt.ylabel("Signed acceleration of the cart (cm/s$^{2}$)")
plt.title("Graph $a$-$m$")
plt.legend()
plt.savefig(FIGURES_DIR / "am.png", dpi=300)
plt.close()

# =========================
# Table 3: linear fit of a_- (negated) against m -- excludes the two zero-acceleration rows
# =========================
fit_minus_mask = is_minus & (mean_a != 0)
m_minus = mass[fit_minus_mask]
a_minus_signed = -mean_a[fit_minus_mask]
u_minus = u_a[fit_minus_mask]

s_minus, u_s_minus, b_minus, u_b_minus, m_minus_x, u_m_minus_x, r2_minus = linear_fit_with_xintercept(
    m_minus, a_minus_signed, u_minus
)

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Slope (\\unit{\\cm\\per\\s\\squared\\per\\g})", "y-intercept (\\unit{\\cm\\per\\s\\squared})", "x-intercept $m_{-}$ (\\unit{\\g})"],
        "Value": [pm(s_minus, u_s_minus), pm(b_minus, u_b_minus), pm(m_minus_x, u_m_minus_x)],
    }),
    caption="Linear regression of $a_{-}$ (sign-flipped) against $m$",
    path=REPORT,
)

x_fit = np.linspace(min(m_minus), max(m_minus), 100)
plot_scatter_fit(
    m_minus, a_minus_signed, u_minus, x_fit, linear_model(x_fit, s_minus, b_minus),
    f'Best-fit line: y = {s_minus:.4g} x - {abs(b_minus):.4g}',
    "Hanging mass $m$ (g)", "Acceleration $a_{-}$ (cm/s$^{2}$)", "Graph $a_{-}$-$m$ with best-fit line",
    "a-withline.png",
)

# =========================
# Table 4: linear fit of a_+ against m -- first 15 (smallest-m) points only
# =========================
m_plus_all = mass[is_plus]
a_plus_all = mean_a[is_plus]
u_plus_all = u_a[is_plus]

m_plus_lin = m_plus_all[:N_FIT_A_PLUS_LINEAR]
a_plus_lin = a_plus_all[:N_FIT_A_PLUS_LINEAR]
u_plus_lin = u_plus_all[:N_FIT_A_PLUS_LINEAR]

s_plus, u_s_plus, b_plus, u_b_plus, m_plus_x, u_m_plus_x, r2_plus = linear_fit_with_xintercept(
    m_plus_lin, a_plus_lin, u_plus_lin
)

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Slope (\\unit{\\cm\\per\\s\\squared\\per\\g})", "y-intercept (\\unit{\\cm\\per\\s\\squared})", "x-intercept $m_{+}$ (\\unit{\\g})"],
        "Value": [pm(s_plus, u_s_plus), pm(b_plus, u_b_plus), pm(m_plus_x, u_m_plus_x)],
    }),
    caption="Linear regression of $a_{+}$ against $m$ (first 15 points only, to avoid the large-$m$ nonlinearity)",
    path=REPORT,
)

x_fit = np.linspace(min(m_plus_lin), max(m_plus_lin), 100)
plot_scatter_fit(
    m_plus_lin, a_plus_lin, u_plus_lin, x_fit, linear_model(x_fit, s_plus, b_plus),
    f'Best-fit line: y = {s_plus:.4g} x - {abs(b_plus):.4g}',
    "Hanging mass $m$ (g)", "Acceleration $a_{+}$ (cm/s$^{2}$)", "Part of graph $a_{+}$-$m$ with its best-fit line",
    "a+withline.png",
)

# =========================
# Table 6: cart mass, friction, and incline angle -- incline-motion method vs. an independent method
# =========================
M_obs = G_CGS / s_minus
u_M_obs = G_CGS * u_s_minus / s_minus**2

Ff_obs_dyne = (m_plus_x - m_minus_x) * G_CGS / 2
u_Ff_obs_dyne = (G_CGS / 2) * np.sqrt(u_m_plus_x**2 + u_m_minus_x**2)
Ff_obs_N = Ff_obs_dyne / DYNE_PER_N
u_Ff_obs_N = u_Ff_obs_dyne / DYNE_PER_N

theta_obs_deg = np.degrees(np.arcsin((m_plus_x + m_minus_x) * s_minus / (2 * G_CGS)))

a_h = np.concatenate([df_horizontal["left_cm_s2"].to_numpy(), df_horizontal["right_cm_s2"].to_numpy()])
a_h_mean, u_a_h = accel_stats(a_h)
Ff_other_N = (M_CART_REF_G / 1000) * abs(a_h_mean) / 100  # g->kg, cm/s^2->m/s^2
u_Ff_other_N = (M_CART_REF_G / 1000) * u_a_h / 100

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Cart mass (\\unit{\\g})", "Friction force (\\unit{\\N})", "Incline angle (\\unit{\\degree})"],
        "Incline-motion method": [pm(M_obs, u_M_obs), pm(Ff_obs_N, u_Ff_obs_N), f"${theta_obs_deg:.3f}$"],
        "Other method": [f"${M_CART_REF_G:.2f}$", pm(Ff_other_N, u_Ff_other_N), f"${theta_other_deg:.3f}$"],
        "Error": [
            percent((M_obs - M_CART_REF_G) * 100 / M_CART_REF_G),
            percent((Ff_obs_N - Ff_other_N) * 100 / Ff_other_N),
            percent((theta_obs_deg - theta_other_deg) * 100 / theta_other_deg),
        ],
    }),
    caption="Comparison of the cart mass, friction force, and incline angle -- incline-motion method vs. an independent method",
    path=REPORT,
)

# =========================
# Tables 11 & 12: nonlinear curve fit of a_pm against m, fixed theta from the trig method
# =========================
def model_a(x, M, F, sign):
    return (G_CGS * x - G_CGS * sin_theta_other * M + sign * F) / (x + M)


popt_minus, pcov_minus = curve_fit(
    lambda x, M, F: model_a(x, M, F, +1), m_minus, a_minus_signed,
    sigma=u_minus, absolute_sigma=True, p0=[M_CART_REF_G, 1500],
)
M_minus_curve, F_minus_curve = popt_minus
u_M_minus_curve, u_F_minus_curve = np.sqrt(np.diag(pcov_minus))

df_to_latex(
    pd.DataFrame({
        "Quantity": ["$M$ (\\unit{\\g})", "$F$ (dyne)"],
        "Value": [pm(M_minus_curve, u_M_minus_curve), pm(F_minus_curve, u_F_minus_curve)],
    }),
    caption="Nonlinear curve fit of $a_{-}$ against $m$",
    path=REPORT,
)

x_fit = np.linspace(min(m_minus), max(m_minus), 100)
plot_scatter_fit(
    m_minus, a_minus_signed, u_minus, x_fit, model_a(x_fit, M_minus_curve, F_minus_curve, +1),
    "Best-fit curve",
    "Hanging mass $m$ (g)", "Acceleration $a_{-}$ (cm/s$^{2}$)", "Graph $a_{-}$-$m$ with best-fit curve",
    "a-withcurve.png",
)

popt_plus, pcov_plus = curve_fit(
    lambda x, M, F: model_a(x, M, F, -1), m_plus_all, a_plus_all,
    sigma=u_plus_all, absolute_sigma=True, p0=[M_CART_REF_G, 900],
)
M_plus_curve, F_plus_curve = popt_plus
u_M_plus_curve, u_F_plus_curve = np.sqrt(np.diag(pcov_plus))

df_to_latex(
    pd.DataFrame({
        "Quantity": ["$M$ (\\unit{\\g})", "$F$ (dyne)"],
        "Value": [pm(M_plus_curve, u_M_plus_curve), pm(F_plus_curve, u_F_plus_curve)],
    }),
    caption="Nonlinear curve fit of $a_{+}$ against $m$",
    path=REPORT,
)

x_fit = np.linspace(min(m_plus_all), max(m_plus_all), 100)
plot_scatter_fit(
    m_plus_all, a_plus_all, u_plus_all, x_fit, model_a(x_fit, M_plus_curve, F_plus_curve, -1),
    "Best-fit curve",
    "Hanging mass $m$ (g)", "Acceleration $a_{+}$ (cm/s$^{2}$)", "Graph $a_{+}$-$m$ with best-fit curve",
    "a+withcurve.png",
)

# =========================
# Table 13: cart mass and friction -- curve-fit (averaged) vs. an independent method
# =========================
M_curve = (M_minus_curve + M_plus_curve) / 2
u_M_curve = 0.5 * np.sqrt(u_M_minus_curve**2 + u_M_plus_curve**2)

Ff_curve_N = (abs(F_minus_curve) + abs(F_plus_curve)) / 2 / DYNE_PER_N
u_Ff_curve_N = 0.5 * np.sqrt(u_F_minus_curve**2 + u_F_plus_curve**2) / DYNE_PER_N

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Cart mass (\\unit{\\g})", "Friction force (\\unit{\\N})"],
        "Curve-fit method": [pm(M_curve, u_M_curve), pm(Ff_curve_N, u_Ff_curve_N)],
        "Other method": [f"${M_CART_REF_G:.2f}$", pm(Ff_other_N, u_Ff_other_N)],
        "Error": [
            percent((M_curve - M_CART_REF_G) * 100 / M_CART_REF_G),
            percent((Ff_curve_N - Ff_other_N) * 100 / Ff_other_N),
        ],
    }),
    caption="Comparison of the cart mass and friction force -- curve-fit method vs. an independent method",
    path=REPORT,
)

# =========================
# Appendix: raw data tables (Tables 14, 15, 16) -- measurement number dropped
# =========================
df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": [f"${v:.2f}$" for v in mass[is_minus]],
        "Acceleration readings $a$ (\\unit{\\cm\\per\\s\\squared})": df_incline.loc[is_minus, "acceleration_cm_s2"],
    }),
    caption="Raw measurements -- incline motion, $a_{-}$ group (3 trials per hanging mass, comma-separated)",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": [f"${v:.2f}$" for v in mass[is_plus]],
        "Acceleration readings $a$ (\\unit{\\cm\\per\\s\\squared})": df_incline.loc[is_plus, "acceleration_cm_s2"],
    }),
    caption="Raw measurements -- incline motion, $a_{+}$ group (3 trials per hanging mass, comma-separated)",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Cart moving left (\\unit{\\cm\\per\\s\\squared})": df_horizontal["left_cm_s2"].map(lambda x: f"${x:.2f}$"),
        "Cart moving right (\\unit{\\cm\\per\\s\\squared})": df_horizontal["right_cm_s2"].map(lambda x: f"${x:.2f}$"),
    }),
    caption="Raw measurements -- cart sliding left/right on horizontal track",
    path=REPORT,
)
