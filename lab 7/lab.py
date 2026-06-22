import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty, linear_model, compute_r_squared

g = 9.8  # m/s^2, standard gravity
L0_cm = 7.50  # cm, spring's natural (unloaded, vertically hung) length
M_g = 11.46   # g, spring mass

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 7" / "figures"
REPORT = BASE_DIR / "lab 7" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

DATA_PATH = BASE_DIR / "lab 7" / "data.xlsx"

df_density = pd.read_excel(DATA_PATH, sheet_name="density")
df_resonance = pd.read_excel(DATA_PATH, sheet_name="resonance")
df_ring = pd.read_excel(DATA_PATH, sheet_name="ring")
df_spring_K = pd.read_excel(DATA_PATH, sheet_name="spring_K")
df_spring_n = pd.read_excel(DATA_PATH, sheet_name="spring_n")


def plot_fit(x, y, y_err, s, b, r_squared, xlabel, ylabel, title, filename):
    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    if y_err is not None:
        plt.errorbar(x, y, yerr=y_err, fmt="o", capsize=3, label="Data Point")
    else:
        plt.scatter(x, y, label="Data Point")
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = linear_model(x_fit, s, b)
    if b >= 0:
        fit_label = f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {r_squared:.4f})'
    else:
        fit_label = f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {r_squared:.4f})'
    plt.plot(x_fit, y_fit, color='black', linestyle='--', label=fit_label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(FIGURES_DIR / filename, dpi=300)
    plt.close()


# =========================
# Part 1: transverse wave on a string -- v = sqrt(T/mu)
# =========================
df_to_latex(
    pd.DataFrame({
        "Mass $m$ (\\unit{\\g})": df_density["mass_g"].map(lambda x: f"${x:.4f}$"),
        "Length $L$ (\\unit{\\cm})": df_density["length_cm"].map(lambda x: f"${x:.2f}$"),
    }),
    path=REPORT,
    caption="Raw mass and length measurements of the string segment",
)

mu_trials = (df_density["mass_g"] / 1000) / (df_density["length_cm"] / 100)
mu_ref, u_mu_ref = typeA_uncertainty(mu_trials.to_numpy())

df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": df_resonance["mass_g"].map(lambda x: f"${x:.2f}$"),
        "Antinode count $n$": df_resonance["n"],
        "Resonance frequency $f_{n}$ (\\unit{\\Hz})": df_resonance["f_Hz"].map(lambda x: f"${x:.1f}$"),
        "Antinode width $\\lambda_{n}/2$ (\\unit{\\cm})": df_resonance["halfwidth_cm"].map(lambda x: f"${x:.1f}$"),
    }),
    path=REPORT,
    caption="Raw standing-wave resonance measurements at different hanging masses",
)

rows1 = []
for mass_g, sub in df_resonance.groupby("mass_g", sort=False):
    v_n = sub["f_Hz"].to_numpy() * 2 * (sub["halfwidth_cm"].to_numpy() / 100)
    v_mean, u_v = typeA_uncertainty(v_n)
    rows1.append({"mass_g": mass_g, "T": mass_g / 1000 * g, "v": v_mean, "u_v": u_v})
part1 = pd.DataFrame(rows1)
part1["v2"] = part1["v"]**2
part1["u_v2"] = 2 * part1["v"] * part1["u_v"]

df_to_latex(
    pd.DataFrame({
        "Rope tension $T=mg$ (\\unit{\\N})": part1["T"].map(lambda x: f"${x:.4f}$"),
        "Wave speed $v=f_{n}\\lambda_{n}$ (\\unit{\\m\\per\\s})": [pm(v, u) for v, u in zip(part1["v"], part1["u_v"])],
    }),
    path=REPORT,
    caption="Wave speed inferred from the resonance measurements, at different rope tensions",
)

popt1, pcov1 = curve_fit(linear_model, part1["T"], part1["v2"], sigma=part1["u_v2"], absolute_sigma=True)
s1, b1 = popt1
u_s1, u_b1 = np.sqrt(np.diag(pcov1))
r2_1 = compute_r_squared(part1["T"], s1, b1, part1["v2"])

plot_fit(part1["T"], part1["v2"], part1["u_v2"], s1, b1, r2_1,
          "Rope tension T (N)", "Squared wave speed v$^{2}$ (m$^{2}$/s$^{2}$)",
          "Graph v$^{2}$-T", "v2T.png")

mu_obs = 1 / s1
u_mu_obs = u_s1 / s1**2
error_mu = (mu_obs - mu_ref) * 100 / mu_ref

df_to_latex(
    pd.DataFrame({
        "Quantity": [
            "Slope $s$ (\\unit{\\m\\squared\\per\\s\\squared\\per\\N})",
            "$\\mu_{\\text{obs}}\\equiv1/s$ (\\unit{\\kg\\per\\m})",
            "$\\mu_{\\text{ref}}$ (\\unit{\\kg\\per\\m})",
            "Error",
        ],
        "Value": [
            pm(s1, u_s1),
            pm(mu_obs, u_mu_obs),
            pm(mu_ref, u_mu_ref),
            percent(error_mu),
        ],
    }),
    path=REPORT,
    caption="Regression of $v^{2}$ against $T$, and the resulting estimate of the string's linear density",
)

# =========================
# Part 2: ring standing wave
# =========================
df_to_latex(
    pd.DataFrame({
        "Antinode count $n$": df_ring["n"],
        "Resonance frequency $f$ (\\unit{\\Hz})": df_ring["f_Hz"].map(lambda x: f"${x:.1f}$"),
    }),
    path=REPORT,
    caption="Raw resonance frequency measurements of the metal ring at different antinode counts",
)

fit_mask = df_ring["n"].isin([3, 5, 7, 9])
n_fit = df_ring.loc[fit_mask, "n"].to_numpy(dtype=float)
f_fit = df_ring.loc[fit_mask, "f_Hz"].to_numpy(dtype=float)

log_n2, log_f2 = np.log10(n_fit), np.log10(f_fit)
popt2, pcov2 = curve_fit(linear_model, log_n2, log_f2)
s2, b2 = popt2
u_s2, u_b2 = np.sqrt(np.diag(pcov2))
r2_2 = compute_r_squared(log_n2, s2, b2, log_f2)

plot_fit(log_n2, log_f2, None, s2, b2, r2_2,
          r"$\log_{10} n$", r"$\log_{10} f$", "Ring standing wave: $\\log f$-$\\log n$",
          "ring_logflogn.png")

error_s2 = (s2 - 2) * 100 / 2
df_to_latex(
    pd.DataFrame({
        "Quantity": ["Slope $s$", "Error"],
        "Value": [pm(s2, u_s2), percent(error_s2)],
    }),
    path=REPORT,
    caption="Regression of $\\log f$ against $\\log n$ for the metal ring, compared with the theoretical exponent of 2",
)

n2_fit = n_fit**2
popt3, pcov3 = curve_fit(linear_model, n2_fit, f_fit)
s3, b3 = popt3
r2_3 = compute_r_squared(n2_fit, s3, b3, f_fit)

plot_fit(n2_fit, f_fit, None, s3, b3, r2_3,
          "$n^{2}$", "Resonance frequency $f$ (Hz)", "Ring standing wave: $f$-$n^{2}$",
          "ring_fn2.png")

rows10 = []
for n_val in [4, 11]:
    f_obs = df_ring.loc[df_ring["n"] == n_val, "f_Hz"].iloc[0]
    f_theo = s3 * n_val**2 + b3
    rows10.append({"n": n_val, "f_obs": f_obs, "f_theo": f_theo,
                    "error": (f_obs - f_theo) * 100 / f_theo})
table10 = pd.DataFrame(rows10)

df_to_latex(
    pd.DataFrame({
        "Antinode count $n$": table10["n"],
        "Observed frequency $f_{\\text{obs}}$ (\\unit{\\Hz})": table10["f_obs"].map(lambda x: f"${x:.1f}$"),
        "Predicted frequency $f_{\\text{theo}}$ (\\unit{\\Hz})": table10["f_theo"].map(lambda x: f"${x:.1f}$"),
        "Error": table10["error"].map(percent),
    }),
    path=REPORT,
    caption="Comparing the observed resonance frequency with the $f$-$n^{2}$ fit's prediction, at antinode counts not used in the fit",
)

# =========================
# Part 3: spring longitudinal wave
# =========================
M_kg = M_g / 1000

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Natural length $L_{0}$ (\\unit{\\cm})", "Spring mass $M$ (\\unit{\\g})"],
        "Value": [f"${L0_cm:.2f}$", f"${M_g:.2f}$"],
    }),
    path=REPORT,
    caption="Basic measurements of the spring used in the longitudinal-wave experiment",
)

df_to_latex(
    pd.DataFrame({
        "Hanging mass $m$ (\\unit{\\g})": df_spring_K["mass_g"].map(lambda x: f"${x:.2f}$"),
        "Upper-loop scale reading $L_{\\text{up}}$ (\\unit{\\cm})": df_spring_K["L_high_cm"].map(lambda x: f"${x:.2f}$"),
        "Lower-loop scale reading $L_{\\text{low}}$ (\\unit{\\cm})": df_spring_K["L_low_cm"].map(lambda x: f"${x:.2f}$"),
    }),
    path=REPORT,
    caption="Raw scale readings used to measure the spring constant",
)

F = df_spring_K["mass_g"].to_numpy() / 1000 * g
x_ext = (df_spring_K["L_high_cm"].to_numpy() - df_spring_K["L_low_cm"].to_numpy() - L0_cm) / 100

df_to_latex(
    pd.DataFrame({
        "Restoring force $F=mg$ (\\unit{\\N})": [f"${v:.4g}$" for v in F],
        "Extension $x=L_{\\text{up}}-L_{\\text{low}}-L_{0}$ (\\unit{\\m})": [f"${v:.4f}$" for v in x_ext],
    }),
    path=REPORT,
    caption="Extension of the spring under different restoring forces",
)

popt4, pcov4 = curve_fit(linear_model, x_ext, F)
K, b4 = popt4
u_K, u_b4 = np.sqrt(np.diag(pcov4))
r2_4 = compute_r_squared(x_ext, K, b4, F)

plot_fit(x_ext, F, None, K, b4, r2_4,
          "Extension $x$ (m)", "Restoring force $F$ (N)", "Graph $F$-$x$",
          "Fx.png")

df_to_latex(
    pd.DataFrame({
        "Quantity": ["Slope $K$ (\\unit{\\N\\per\\m})", "y-intercept (\\unit{\\N})"],
        "Value": [pm(K, u_K), pm(b4, u_b4)],
    }),
    path=REPORT,
    caption="Regression of $F$ against $x$; by Hooke's law the slope is the spring constant $K$",
)

df_to_latex(
    pd.DataFrame({
        "Antinode count $n$": df_spring_n["n"],
        "Resonance frequency $f$ (\\unit{\\Hz})": df_spring_n["f_Hz"].map(lambda x: f"${x:.1f}$"),
    }),
    path=REPORT,
    caption="Raw resonance frequency measurements of the spring's longitudinal standing wave",
)

n5 = df_spring_n["n"].to_numpy(dtype=float)
f5 = df_spring_n["f_Hz"].to_numpy(dtype=float)

popt5, pcov5 = curve_fit(linear_model, n5, f5)
s5, b5 = popt5
u_s5, u_b5 = np.sqrt(np.diag(pcov5))
r2_5 = compute_r_squared(n5, s5, b5, f5)

plot_fit(n5, f5, None, s5, b5, r2_5,
          "Antinode count $n$", "Resonance frequency $f$ (Hz)", "Spring longitudinal wave: $f$-$n$",
          "spring_fn.png")

s_theo = np.sqrt(K / M_kg) / 2
u_s_theo = s_theo / 2 / K * u_K
error_s5 = (s5 - s_theo) * 100 / s_theo

df_to_latex(
    pd.DataFrame({
        "Quantity": [
            "Slope $s$ (\\unit{\\Hz})",
            "Theoretical slope $s_{\\text{theo}}=\\sqrt{K/M}/2$ (\\unit{\\Hz})",
            "Error",
        ],
        "Value": [pm(s5, u_s5), pm(s_theo, u_s_theo), percent(error_s5)],
    }),
    path=REPORT,
    caption="Regression of $f$ against $n$ for the spring, compared with the theoretical slope $\\sqrt{K/M}/2$",
)

log_n5, log_f5 = np.log10(n5), np.log10(f5)
popt6, pcov6 = curve_fit(linear_model, log_n5, log_f5)
s6, b6 = popt6
u_s6, u_b6 = np.sqrt(np.diag(pcov6))
r2_6 = compute_r_squared(log_n5, s6, b6, log_f5)

plot_fit(log_n5, log_f5, None, s6, b6, r2_6,
          r"$\log_{10} n$", r"$\log_{10} f$", "Spring longitudinal wave: $\\log f$-$\\log n$",
          "spring_logflogn.png")

error_s6 = (s6 - 1) * 100 / 1
b_theo = np.log10(s_theo)
u_b_theo = u_s_theo / (s_theo * np.log(10))
error_b6 = (b6 - b_theo) * 100 / b_theo

df_to_latex(
    pd.DataFrame({
        "Quantity": [
            "Slope $s$",
            "Error",
            "y-intercept $b$",
            "Theoretical y-intercept $b_{\\text{theo}}=\\log s_{\\text{theo}}$",
            "Error",
        ],
        "Value": [
            pm(s6, u_s6),
            percent(error_s6),
            pm(b6, u_b6),
            pm(b_theo, u_b_theo),
            percent(error_b6),
        ],
    }),
    path=REPORT,
    caption="Regression of $\\log f$ against $\\log n$ for the spring, compared with the theoretical exponent of 1",
)
