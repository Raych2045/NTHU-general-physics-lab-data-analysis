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
g = 9.8  # m/s^2, standard gravity

# Parts 1-3 share the same disk and tower pulley
r_pulley = 0.0125     # m, tower-pulley radius (2r = 25.00 mm), fixed in parts 1-3
M_disk = 1.4453       # kg, disk mass (parts 1-3)
R_disk = 0.1142       # m, disk radius (2R = 22.84 cm)
L_disk = 0.0254       # m, disk thickness (used in Part 2's reference moment of inertia)

M_disk_g = M_disk * 1000
R_disk_diam_cm = R_disk * 2 * 100
L_disk_mm = L_disk * 1000

I_hori_ref = M_disk * R_disk**2 / 2
I_verti_ref = M_disk * R_disk**2 / 4 + M_disk * L_disk**2 / 12
mass_ref_kg = 1.4757  # disk + adapter mass, compared to the I_off-d^2 fit's slope

T_TYPEB_MS = 1 / np.sqrt(12)  # ms, Type B uncertainty of each T1/T2 timer reading (1 ms resolution)

# Part 4: torsion pendulum, different (reference) disk used to calibrate kappa
M_ref4 = 0.6606          # kg, calibration disk mass
R_ref4 = 0.12030 / 2     # m, calibration disk radius (2R = 120.30 mm)
I_ref4 = M_ref4 * R_ref4**2 / 2

M_ref4_g = M_ref4 * 1000
R_ref4_diam_mm = R_ref4 * 2 * 1000

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 5" / "figures"
REPORT = BASE_DIR / "lab 5" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

DATA_PATH = BASE_DIR / "lab 5" / "data.xlsx"

# =========================
# load data
# =========================
df_ts = pd.read_excel(DATA_PATH, sheet_name="part1_3_timestamps")
df_cond = pd.read_excel(DATA_PATH, sheet_name="part1_3_conditions").set_index("group")
df4_raw = pd.read_excel(DATA_PATH, sheet_name="part4_raw")

def _group_label(group):
    if group == "hori":
        return "Horizontal disk (about its central axis)"
    if group == "verti":
        return "Vertical disk (about its diameter axis)"
    if group == "platform":
        return "Rotating platform + counterweight only (no disk)"
    d_cm = df_cond.loc[group, "d_cm"]
    return f"Off-axis disk at $d = {d_cm:.2f}$ cm"


GROUP_LABELS = {group: _group_label(group) for group in df_cond.index}


# =========================
# Part 1-3 helper: weighted f-t linear fit -> angular acceleration
# =========================
def ft_fit(run_name):
    chain = df_ts[run_name].dropna().to_numpy(dtype=float)
    T1, T2 = chain[:-1], chain[1:]
    delta = T2 - T1
    t = (T1 + T2) / 2000.0           # midpoint time, s (T1,T2 in ms)
    f = 100.0 / delta                # Hz, since f = 0.1 / ((T2-T1)/1000 s)
    u_f = 100.0 * np.sqrt(1 / 6) / delta**2

    popt, pcov = curve_fit(linear_model, t, f, sigma=u_f, absolute_sigma=True, p0=[0.0, 0.0])
    s, b = popt
    u_s, u_b = np.sqrt(np.diag(pcov))
    r_squared = compute_r_squared(t, s, b, f)
    return {"t": t, "f": f, "u_f": u_f, "s": s, "u_s": u_s, "b": b, "u_b": u_b, "r_squared": r_squared}


def plot_ft(fit, group, run_kind, label):
    s, b, r_squared = fit["s"], fit["b"], fit["r_squared"]
    t_fit = np.linspace(min(fit["t"]), max(fit["t"]), 100)
    f_fit = linear_model(t_fit, s, b)

    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.errorbar(fit["t"], fit["f"], yerr=fit["u_f"], fmt="o", capsize=3, label="Data Point")
    if b >= 0:
        fit_label = f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {r_squared:.4f})'
    else:
        fit_label = f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {r_squared:.4f})'
    plt.plot(t_fit, f_fit, color='black', linestyle='--', label=fit_label)
    plt.xlabel(r"Midpoint time $t = (T_1+T_2)/2$ (s)")
    plt.ylabel(r"Frequency $f = 0.1/(T_2-T_1)$ (Hz)")
    plt.title(f"{label} -- {run_kind} run")
    plt.legend()
    plt.savefig(FIGURES_DIR / f"{group}_{run_kind}.png", dpi=300)
    plt.close()


def moment_of_inertia(m_kg, alpha, u_alpha, alpha_f, u_alpha_f):
    I = m_kg * r_pulley * (g - r_pulley * alpha) / (alpha - alpha_f)
    front = m_kg * r_pulley / (alpha - alpha_f)**2
    mid = r_pulley * alpha_f - g
    back = g - r_pulley * alpha
    u_I = front * np.sqrt((mid * u_alpha)**2 + (back * u_alpha_f)**2)
    return I, u_I


results = {}
for group in GROUP_LABELS:
    loaded = ft_fit(f"{group}_loaded")
    friction = ft_fit(f"{group}_friction")
    plot_ft(loaded, group, "loaded", GROUP_LABELS[group])
    plot_ft(friction, group, "friction", GROUP_LABELS[group])

    alpha, u_alpha = 2 * np.pi * loaded["s"], 2 * np.pi * loaded["u_s"]
    alpha_f, u_alpha_f = 2 * np.pi * friction["s"], 2 * np.pi * friction["u_s"]
    m_kg = df_cond.loc[group, "mass_g"] / 1000
    I, u_I = moment_of_inertia(m_kg, alpha, u_alpha, alpha_f, u_alpha_f)

    results[group] = {
        "alpha": alpha, "u_alpha": u_alpha,
        "alpha_f": alpha_f, "u_alpha_f": u_alpha_f,
        "I": I, "u_I": u_I,
    }

# =========================
# Part 1: horizontal disk
# =========================
r1 = results["hori"]
error1 = (r1["I"] - I_hori_ref) * 100 / I_hori_ref
table1 = pd.DataFrame({
    "Quantity": [
        "$\\alpha$ (\\unit{\\radian\\per\\s\\squared})",
        "$\\alpha_{f}$ (\\unit{\\radian\\per\\s\\squared})",
        "$I_{\\text{hori,obs}}$ (\\unit{\\kg\\m\\squared})",
        "Disk mass $M$ (\\unit{\\g})",
        "Disk diameter $2R$ (\\unit{\\cm})",
        "$I_{\\text{hori,ref}}=MR^{2}/2$ (\\unit{\\kg\\m\\squared})",
        "Error",
    ],
    "Value": [
        pm(r1["alpha"], r1["u_alpha"]),
        pm(r1["alpha_f"], r1["u_alpha_f"]),
        pm(r1["I"], r1["u_I"]),
        f"${M_disk_g:.1f}$",
        f"${R_disk_diam_cm:.2f}$",
        f"${I_hori_ref:.6f}$",
        percent(error1),
    ],
})
df_to_latex(table1, path=REPORT,
            caption="Moment of inertia of the disk about its central axis")

# =========================
# Part 2: vertical disk
# =========================
r2 = results["verti"]
error2 = (r2["I"] - I_verti_ref) * 100 / I_verti_ref
table2 = pd.DataFrame({
    "Quantity": [
        "$\\alpha$ (\\unit{\\radian\\per\\s\\squared})",
        "$\\alpha_{f}$ (\\unit{\\radian\\per\\s\\squared})",
        "$I_{\\text{verti,obs}}$ (\\unit{\\kg\\m\\squared})",
        "Disk thickness $L$ (\\unit{\\mm})",
        "$I_{\\text{verti,ref}}=MR^{2}/4+ML^{2}/12$ (\\unit{\\kg\\m\\squared})",
        "Error",
    ],
    "Value": [
        pm(r2["alpha"], r2["u_alpha"]),
        pm(r2["alpha_f"], r2["u_alpha_f"]),
        pm(r2["I"], r2["u_I"]),
        f"${L_disk_mm:.2f}$",
        f"${I_verti_ref:.6f}$",
        percent(error2),
    ],
})
df_to_latex(table2, path=REPORT,
            caption="Moment of inertia of the disk about its diameter axis")

# =========================
# Part 3: off-axis disk (parallel axis theorem) -- one small table per d / platform
# =========================
for group in ["d10", "d11", "d12", "d13", "d14"]:
    r = results[group]
    d_cm = df_cond.loc[group, "d_cm"]
    table_d = pd.DataFrame({
        "Quantity": [
            "$\\alpha$ (\\unit{\\radian\\per\\s\\squared})",
            "$\\alpha_{f}$ (\\unit{\\radian\\per\\s\\squared})",
            "$I_{\\text{off}}$ (\\unit{\\kg\\m\\squared})",
        ],
        "Value": [
            pm(r["alpha"], r["u_alpha"]),
            pm(r["alpha_f"], r["u_alpha_f"]),
            pm(r["I"], r["u_I"]),
        ],
    })
    df_to_latex(table_d, path=REPORT,
                caption=f"Moment of inertia of the whole system with the disk offset by $d=\\qty{{{d_cm:.2f}}}{{\\cm}}$")

r_plat = results["platform"]
table_plat = pd.DataFrame({
    "Quantity": [
        "$\\alpha$ (\\unit{\\radian\\per\\s\\squared})",
        "$\\alpha_{f}$ (\\unit{\\radian\\per\\s\\squared})",
        "$I_{\\text{plat}}$ (\\unit{\\kg\\m\\squared})",
    ],
    "Value": [
        pm(r_plat["alpha"], r_plat["u_alpha"]),
        pm(r_plat["alpha_f"], r_plat["u_alpha_f"]),
        pm(r_plat["I"], r_plat["u_I"]),
    ],
})
df_to_latex(table_plat, path=REPORT,
            caption="Moment of inertia of the rotating platform and counterweight alone (no disk)")

# Summary table: d, I_off
groups_d = ["d10", "d11", "d12", "d13", "d14"]
d_cm_arr = df_cond.loc[groups_d, "d_cm"].to_numpy(dtype=float)
I_off_arr = np.array([results[g]["I"] for g in groups_d])
u_I_off_arr = np.array([results[g]["u_I"] for g in groups_d])

table12 = pd.DataFrame({
    "Offset $d$ (\\unit{\\cm})": [f"${d:.2f}$" for d in d_cm_arr],
    "$I_{\\text{off}}$ (\\unit{\\kg\\m\\squared})": [pm(i, u) for i, u in zip(I_off_arr, u_I_off_arr)],
})
df_to_latex(table12, path=REPORT,
            caption="Moment of inertia of the whole system $I_{\\text{off}}$ at five different disk offsets $d$")

# Parallel axis theorem regression: I_off vs d^2
d2_arr = (d_cm_arr / 100)**2  # m^2
popt3, pcov3 = curve_fit(linear_model, d2_arr, I_off_arr, sigma=u_I_off_arr, absolute_sigma=True)
s3, b3 = popt3
u_s3, u_b3 = np.sqrt(np.diag(pcov3))
r_squared3 = compute_r_squared(d2_arr, s3, b3, I_off_arr)

d2_fit = np.linspace(min(d2_arr), max(d2_arr), 100)
I_fit = linear_model(d2_fit, s3, b3)

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(d2_arr, I_off_arr, yerr=u_I_off_arr, fmt="o", capsize=3, label="Data Point")
if b3 >= 0:
    fit_label3 = f'Best-fit line: y = {s3:.4g} x + {b3:.4g} (R² = {r_squared3:.4f})'
else:
    fit_label3 = f'Best-fit line: y = {s3:.4g} x - {abs(b3):.4g} (R² = {r_squared3:.4f})'
plt.plot(d2_fit, I_fit, color='black', linestyle='--', label=fit_label3)
plt.xlabel(r"$d^{2}$ (m$^{2}$)")
plt.ylabel(r"$I_{\text{off}}$ (kg m$^{2}$)")
plt.title("Parallel axis theorem: $I_{\\text{off}}$ versus $d^{2}$")
plt.legend()
plt.savefig(FIGURES_DIR / "parallel_axis.png", dpi=300)
plt.close()

mass_ref_g = mass_ref_kg * 1000
slope_obs_g = s3 * 1000
u_slope_obs_g = u_s3 * 1000
error_slope = (s3 - mass_ref_kg) * 100 / mass_ref_kg

I_plat, u_I_plat = r_plat["I"], r_plat["u_I"]
b_minus_Iplat = b3 - I_plat
u_b_minus_Iplat = np.sqrt(u_b3**2 + u_I_plat**2)
error_intercept = (b_minus_Iplat - I_hori_ref) * 100 / I_hori_ref

table13 = pd.DataFrame({
    "Quantity": [
        "Slope $s$ (\\unit{\\g})",
        "Disk + adapter mass $M_{\\text{ref}}$ (\\unit{\\g})",
        "Error",
        "y-intercept $b$ (\\unit{\\kg\\m\\squared})",
        "$b-I_{\\text{plat}}$ (\\unit{\\kg\\m\\squared})",
        "$I_{\\text{hori,ref}}$ (\\unit{\\kg\\m\\squared})",
        "Error",
    ],
    "Value": [
        pm(slope_obs_g, u_slope_obs_g),
        f"${mass_ref_g:.1f}$",
        percent(error_slope),
        pm(b3, u_b3),
        pm(b_minus_Iplat, u_b_minus_Iplat),
        f"${I_hori_ref:.6f}$",
        percent(error_intercept),
    ],
})
df_to_latex(table13, path=REPORT,
            caption="Parallel axis theorem -- comparing the $I_{\\text{off}}$ - $d^{2}$ fit's slope and intercept with theory")

# =========================
# Part 4: torsion pendulum (black box)
# =========================
def thirty_period_stats(config):
    sub = df4_raw[df4_raw["config"] == config]
    thirty_T = (sub["T2_s"] - sub["T1_s"]).to_numpy(dtype=float)
    mean_30T, u_30T = typeA_uncertainty(thirty_T)
    T, u_T = mean_30T / 30, u_30T / 30
    return mean_30T, u_30T, T, u_T


thirtyT_disk, u_thirtyT_disk, T_disk, u_T_disk = thirty_period_stats("disk")
kappa = 4 * np.pi**2 * I_ref4 / T_disk**2
u_kappa = 8 * np.pi**2 * I_ref4 * u_T_disk / T_disk**3

table14 = pd.DataFrame({
    "Quantity": [
        "$30T$ (\\unit{\\s})",
        "Period $T$ (\\unit{\\s})",
        "Calibration disk mass $M$ (\\unit{\\g})",
        "Calibration disk diameter $2R$ (\\unit{\\mm})",
        "Calibration disk's $I=M R^{2}/2$ (\\unit{\\kg\\m\\squared})",
        "$\\kappa$ (\\unit{\\N\\m\\per\\radian})",
    ],
    "Value": [
        pm(thirtyT_disk, u_thirtyT_disk),
        pm(T_disk, u_T_disk),
        f"${M_ref4_g:.1f}$",
        f"${R_ref4_diam_mm:.2f}$",
        f"${I_ref4:.6f}$",
        pm(kappa, u_kappa),
    ],
})
df_to_latex(table14, path=REPORT,
            caption="Torsion constant $\\kappa$ of the suspension wire, calibrated with a disk of known moment of inertia")

HOLE_LABELS = {"hole1": "I_{xx}", "hole2": "I_{yy}", "hole3": "I_{zz}"}
for hole, sym in HOLE_LABELS.items():
    thirtyT_h, u_thirtyT_h, T_h, u_T_h = thirty_period_stats(hole)
    I_h = kappa * (T_h / (2 * np.pi))**2
    u_I_h = (1 / (4 * np.pi**2)) * np.sqrt(T_h**4 * u_kappa**2 + 4 * kappa**2 * T_h**2 * u_T_h**2)

    table_h = pd.DataFrame({
        "Quantity": [
            "$30T$ (\\unit{\\s})",
            "Period $T$ (\\unit{\\s})",
            f"${sym}$ (\\unit{{\\kg\\m\\squared}})",
        ],
        "Value": [
            pm(thirtyT_h, u_thirtyT_h),
            pm(T_h, u_T_h),
            pm(I_h, u_I_h),
        ],
    })
    df_to_latex(table_h, path=REPORT,
                caption=f"Moment of inertia of the black box about the axis through ${sym}$'s suspension hole")

# =========================
# Raw data (appendix)
# =========================
def raw_chain_table(run_name, label):
    chain = df_ts[run_name].dropna().to_numpy(dtype=float)
    values = [f"{v:.0f}" for v in chain]
    lines = [", ".join(values[i:i + 10]) for i in range(0, len(values), 10)]
    cell = "\\makecell{" + " \\\\ ".join(lines) + "}"
    table = pd.DataFrame({
        "Stopwatch readings $T$ (\\unit{\\ms})": [cell],
    })
    df_to_latex(table, path=REPORT, caption=label)


for group, label in GROUP_LABELS.items():
    raw_chain_table(f"{group}_loaded", f"Raw stopwatch readings -- {label}, loaded run")
    raw_chain_table(f"{group}_friction", f"Raw stopwatch readings -- {label}, friction (unloaded) run")

raw4 = df4_raw.copy()
raw4["30T_s"] = raw4["T2_s"] - raw4["T1_s"]
table_raw4 = pd.DataFrame({
    "Configuration": raw4["config"],
    "$T_1$ (\\unit{\\s})": raw4["T1_s"].map(lambda x: f"${x:.2f}$"),
    "$T_2$ (\\unit{\\s})": raw4["T2_s"].map(lambda x: f"${x:.2f}$"),
    "$30T=T_2-T_1$ (\\unit{\\s})": raw4["30T_s"].map(lambda x: f"${x:.2f}$"),
})
df_to_latex(table_raw4, path=REPORT,
            caption="Raw timing measurements for the calibration disk and the three black-box suspension holes")
