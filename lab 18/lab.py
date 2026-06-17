import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, percent, df_to_latex, typeA_uncertainty

BASE_DIR    = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 18" / "figures"
REPORT      = BASE_DIR / "lab 18" / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")

# =========================
# Constants
# =========================
n_ref          = 1.49
lambda_ref_nm  = 532
lambda_ref_m   = lambda_ref_nm * 1e-9

def rel_err(obs, ref):
    return (obs - ref) / ref * 100.0

def extract(s):
    return np.array([float(x.strip()) for x in str(s).split(";")])

def fmt_trig(val):
    """Format a trig value, snapping to 0 or 1 when within floating-point noise."""
    if abs(val) < 1e-10:
        return "$0$"
    if abs(val - 1.0) < 1e-10:
        return "$1$"
    return f"${val:.4g}$"

# =========================
# Load data
# =========================
xl   = pd.ExcelFile(BASE_DIR / "lab 18" / "data.xlsx")
df_A = xl.parse("A")
df_B = xl.parse("B")
df_C = xl.parse("C")
df_m2 = xl.parse("malus_two")
df_m3 = xl.parse("malus_three")
df_ss = xl.parse("single_slit")
df_ds = xl.parse("double_slit")

# =========================
# Sub-experiment A: acrylic plate
# n = sqrt(d^2 + r^2) / r,  where r = R/2 = (2R)/4
# Only Type-A uncertainty on d; 2R is taken as exact.
# =========================
results_A = []
for desc, grp in df_A.groupby("description"):
    R2_note = grp["R2_note"].iloc[0]
    d_vals  = extract(grp["d_mm"].iloc[0])
    R2      = grp["2R_mm"].iloc[0]

    d_mean, uA_d = typeA_uncertainty(d_vals)
    u_d  = np.sqrt(uA_d**2 + (0.02 / np.sqrt(12))**2)   # type B: 0.02 mm

    r    = R2 / 4                                          # r = R/2
    n    = np.sqrt(d_mean**2 + r**2) / r
    dn_dd = d_mean / (r * np.sqrt(d_mean**2 + r**2))
    u_n  = abs(dn_dd) * u_d

    results_A.append(dict(
        desc=desc, R2_note=R2_note,
        d_vals=d_vals, d_mean=d_mean, u_d=u_d,
        R2=R2, n=n, u_n=u_n,
        err=rel_err(n, n_ref),
    ))

# Table 1
df_to_latex(
    pd.DataFrame({
        "描述": [
            "厚度$d$ (\\unit{\\mm})",
            "$d$的觀測結果 (\\unit{\\mm})",
            "直徑$2R$ (\\unit{\\mm})",
            "$n_{\\text{obs}}$",
            "Error",
        ],
        **{r["desc"]: [
            ", ".join(f"{v:.2f}" for v in r["d_vals"]),
            pm(r["d_mean"], r["u_d"]),
            f"${r['R2']:.2f}$ （{r['R2_note']}）",
            pm(r["n"], r["u_n"]),
            percent(r["err"]),
        ] for r in results_A},
    }),
    caption="壓克力板子實驗的數據，未考慮$2R$的不確定度。",
    path=REPORT,
)

# =========================
# Sub-experiment B: acrylic rod
# n = (D_acry - l) / (D_air - l),  all quantities in metres
# =========================
D_vals   = df_B["D_acry_m"].values
da_vals  = df_B["D_air_m"].values
l_mm_v   = df_B["l_mm"].values

D_mean,  uA_D  = typeA_uncertainty(D_vals)
da_mean, uA_da = typeA_uncertainty(da_vals)
lmm_mean, uA_lmm = typeA_uncertainty(l_mm_v)

u_D  = np.sqrt(uA_D**2  + (0.001 / np.sqrt(12))**2)   # type B: 1 mm
u_da = np.sqrt(uA_da**2 + (0.001 / np.sqrt(12))**2)
u_lmm = np.sqrt(uA_lmm**2 + (0.02 / np.sqrt(12))**2)  # type B: 0.02 mm

l_mean = lmm_mean / 1000   # mm → m
u_l    = u_lmm / 1000

n_B    = (D_mean - l_mean) / (da_mean - l_mean)
dn_dD  =  1 / (da_mean - l_mean)
dn_dda = -(D_mean  - l_mean) / (da_mean - l_mean)**2
dn_dl  =  (D_mean  - da_mean) / (da_mean - l_mean)**2
u_n_B  = np.sqrt((dn_dD * u_D)**2 + (dn_dda * u_da)**2 + (dn_dl * u_l)**2)

# Table 2 (raw data)
df_to_latex(
    pd.DataFrame({
        "壓克力柱$D_{\\text{acry}}$ (\\unit{\\m})": df_B["D_acry_m"].map(lambda x: f"${x:.3f}$"),
        "空氣柱$D_{\\text{air}}$ (\\unit{\\m})":    df_B["D_air_m"].map(lambda x: f"${x:.3f}$"),
        "測距儀$l$ (\\unit{\\mm})":                 df_B["l_mm"].map(lambda x: f"${x:.2f}$"),
    }),
    caption="壓克力棒子實驗的數據",
    path=REPORT,
)

# Table 3 (results)
df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "$D_{\\text{acry}}$的觀測結果 (\\unit{\\m})",
            "$D_{\\text{air}}$的觀測結果 (\\unit{\\m})",
            "$l$的觀測結果 (\\unit{\\mm})",
            "$n_{\\text{obs}}=\\dfrac{D_{\\text{acry}}-l}{D_{\\text{air}}-l}$",
            "Error",
        ],
        "Value": [
            pm(D_mean, u_D),
            pm(da_mean, u_da),
            pm(lmm_mean, u_lmm),
            pm(n_B, u_n_B),
            percent(rel_err(n_B, n_ref)),
        ],
    }),
    caption="續上表。這裡未考慮$D_{\\text{acry}}$和$D_{\\text{air}}$的B類不確定度。",
    path=REPORT,
)

# =========================
# Two-polariser Malus' law
# =========================
theta_2    = df_m2["theta_deg"].values.astype(float)
Iavg_2     = []
u_I_2      = []
raw_str_2  = []

for _, row in df_m2.iterrows():
    vals = extract(row["Iobs"])
    m, u = typeA_uncertainty(vals)
    Iavg_2.append(m)
    u_I_2.append(u)
    raw_str_2.append(", ".join(f"{v:.3f}" for v in vals))

I0       = Iavg_2[0]
cos2_arr = np.cos(np.deg2rad(theta_2))**2
Itheo_2  = I0 * cos2_arr
theta_fit = np.linspace(min(theta_2), max(theta_2), 100)
Itheo_fit = I0 * np.cos(np.deg2rad(theta_fit))**2

# Plot two.png
plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(theta_2, Iavg_2, yerr=u_I_2, fmt="o", capsize=3, label="Data Point")
plt.plot(theta_fit, Itheo_fit, color='black', linestyle='--',
         label="Prediction from Malus' Law")
plt.xlabel("Angle (degree)")
plt.ylabel("Illuminance (lx)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.savefig(FIGURES_DIR / "two.png", dpi=300)
plt.close()

# Table 4
def fmt_Itheo(v):
    return "$0$" if abs(v) < 1e-10 else f"${v:.3f}$"

df_to_latex(
    pd.DataFrame({
        "$\\theta$ (\\unit{\\degree})":                             [f"${int(t)}$" for t in theta_2],
        "$\\cos^{2}\\theta$":                                       [fmt_trig(c) for c in cos2_arr],
        "Raw $I_{\\text{obs}}$ (\\unit{\\lux})":                   raw_str_2,
        "$I_{\\text{obs}}$ (\\unit{\\lux})":                                     [pm(m, u) for m, u in zip(Iavg_2, u_I_2)],
        "$I_{\\text{theo}}:=I_{0}\\cos^{2}\\theta$ (\\unit{\\lux})": [fmt_Itheo(v) for v in Itheo_2],
    }),
    caption="使用兩片偏振片觀測透光強度，其中$I_{0}$取$\\theta=0$時的$I_{\\text{obs}}$。",
    path=REPORT,
)

# Table 5 (back-computed angle, theta = 30, 45, 60, 90 only)
mask       = theta_2 != 0
th_non0    = theta_2[mask]
I_non0     = np.array(Iavg_2)[mask]
th_obs_arr = np.degrees(np.arccos(np.sqrt(np.clip(I_non0 / I0, 0, 1))))

df_to_latex(
    pd.DataFrame({
        " ": ["$\\theta_{\\text{obs}}:=\\arccos\\sqrt{I/I_{0}}$ (\\unit{\\degree})"],
        **{f"${int(t)}$": [f"${v:.1f}$"] for t, v in zip(th_non0, th_obs_arr)},
    }),
    caption="利用光強觀測值反推穿射軸夾角",
    path=REPORT,
)

# =========================
# Three-polariser Malus' law
# =========================
theta_3    = df_m3["theta_deg"].values.astype(float)
Iavg_3     = []
u_I_3      = []
raw_str_3  = []

for _, row in df_m3.iterrows():
    vals = extract(row["Iobs"])
    m, u = typeA_uncertainty(vals)
    Iavg_3.append(m)
    u_I_3.append(u)
    raw_str_3.append(", ".join(f"{v:.3f}" for v in vals))

sin2_2th_3 = np.sin(2 * np.deg2rad(theta_3))**2

# Plot three.png
plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
plt.errorbar(theta_3, Iavg_3, yerr=u_I_3, fmt="o", capsize=3, label="Data Point")
plt.xlabel("Angle (degree)")
plt.ylabel("Illuminance (lx)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.savefig(FIGURES_DIR / "three.png", dpi=300)
plt.close()

# Table 6
df_to_latex(
    pd.DataFrame({
        "$\\theta$ (\\unit{\\degree})":          [f"${int(t)}$" for t in theta_3],
        "$\\sin^{2}2\\theta$":                   [fmt_trig(v) for v in sin2_2th_3],
        "Raw $I_{\\text{obs}}$ (\\unit{\\lux})": raw_str_3,
        "$I_{\\text{obs}}$ (\\unit{\\lux})":                   [pm(m, u) for m, u in zip(Iavg_3, u_I_3)],
    }),
    caption="使用三片偏振片觀測透光強度",
    path=REPORT,
)

# =========================
# Sub-experiment C: Brewster's angle
# n = tan(theta_B)
# =========================
tB_vals  = df_C["theta_B_deg"].values
t_mean, uA_t = typeA_uncertainty(tB_vals)
u_t = np.sqrt(uA_t**2 + (1.0 / np.sqrt(12))**2)   # type B: 1°

t_rad = np.deg2rad(t_mean)
u_rad = np.deg2rad(u_t)
n_C   = np.tan(t_rad)
u_n_C = abs(1 / np.cos(t_rad)**2) * u_rad

# Table 7
df_to_latex(
    pd.DataFrame({
        "Parameter": [
            "$\\theta_{B}$ (\\unit{\\degree})",
            "$\\theta_{B}$的測量結果 (\\unit{\\degree})",
            "$n_{\\text{obs}}=\\tan\\theta_{B}$",
            "Error",
        ],
        "Value": [
            ", ".join(f"{v:.1f}" for v in tB_vals),
            pm(t_mean, u_t),
            pm(n_C, u_n_C),
            percent(rel_err(n_C, n_ref)),
        ],
    }),
    caption="Brewster's angle的觀測結果",
    path=REPORT,
)

# =========================
# Single-slit diffraction
# lambda_obs = a * W / (2 * L),  a and W in m, L in m
# =========================
lam_obs_list = []
err_ss_list  = []

for _, row in df_ss.iterrows():
    a_m   = row["a_mm"] * 1e-3
    W_m   = row["W_mm"] * 1e-3
    L_m   = row["L_cm"] * 1e-2
    lam   = a_m * W_m / (2 * L_m) * 1e9   # → nm
    lam_obs_list.append(lam)
    err_ss_list.append(rel_err(lam, lambda_ref_nm))

L_ss = df_ss["L_cm"].iloc[0]

df_to_latex(
    pd.DataFrame({
        "$a$ (\\unit{\\mm})":                             df_ss["a_mm"].map(lambda x: f"${x:.2f}$"),
        "$W$ (\\unit{\\mm})":                             df_ss["W_mm"].map(lambda x: f"${x:.2f}$"),
        "$\\lambda_{\\text{obs}}=aW/2L$ (\\unit{\\nm})": [f"${v:.0f}$" for v in lam_obs_list],
        "Error":                                          [percent(e) for e in err_ss_list],
    }),
    caption=(
        f"單狹縫繞射子實驗的數據，其中$L= \\qty{{{L_ss:.2f}}}{{\\cm}}$。"
        f"$\\lambda_{{\\text{{ref}}}}$取 \\qty{{{lambda_ref_nm}}}{{\\nm}}。"
    ),
    path=REPORT,
)

# =========================
# Double-slit diffraction
# a_obs = 2 * L * lambda_ref / W
# d_obs = L * lambda_ref / delta_y,  delta_y = span_mm / N_fringes
# =========================
a_obs_list   = []
d_obs_list   = []
err_a_list   = []
err_d_list   = []
dy_str_list  = []

for _, row in df_ds.iterrows():
    L_m      = row["L_cm"] * 1e-2
    W_m      = row["W_mm"] * 1e-3
    span_mm  = float(row["span_mm"])
    N        = int(row["N_fringes"])
    a_ref_mm = float(row["a_ref_mm"])
    d_ref_mm = float(row["d_ref_mm"])

    dy_mm  = span_mm / N
    dy_m   = dy_mm * 1e-3

    a_obs_mm = 2 * L_m * lambda_ref_m / W_m * 1e3
    d_obs_mm = L_m * lambda_ref_m / dy_m * 1e3

    a_obs_list.append(a_obs_mm)
    d_obs_list.append(d_obs_mm)
    err_a_list.append(rel_err(a_obs_mm, a_ref_mm))
    err_d_list.append(rel_err(d_obs_mm, d_ref_mm))
    dy_str_list.append(f"${dy_mm:.3f}$ (={span_mm:.2f}/{N})")

L_ds = df_ds["L_cm"].iloc[0]

# Table 10
df_to_latex(
    pd.DataFrame({
        "代號":                                df_ds["label"].map(str),
        "$W$ (\\unit{\\mm})":                 df_ds["W_mm"].map(lambda x: f"${x:.2f}$"),
        "$a_{\\text{obs}}$ (\\unit{\\mm})":   [f"${v:.3f}$" for v in a_obs_list],
        "$a_{\\text{ref}}$ (\\unit{\\mm})":   df_ds["a_ref_mm"].map(lambda x: f"${x:.2f}$"),
        "Error":                         [percent(e) for e in err_a_list],
        "$\\Delta y$ (\\unit{\\mm})":          dy_str_list,
        "$d_{\\text{obs}}$ (\\unit{\\mm})":   [f"${v:.3f}$" for v in d_obs_list],
        "$d_{\\text{ref}}$ (\\unit{\\mm})":   df_ds["d_ref_mm"].map(lambda x: f"${x:.2f}$"),
        "Error":                         [percent(e) for e in err_d_list],
    }),
    caption=(
        f"雙狹縫繞射子實驗的數據，其中$L= \\qty{{{L_ds:.2f}}}{{\\cm}}$，"
        f"$\\lambda$取$\\lambda_{{\\text{{ref}}}}$"
    ),
    path=REPORT,
)

# Table 11 (back-computed tilt angle phi = arccos(d_obs / d_ref))
phi_list = [
    np.degrees(np.arccos(np.clip(d_obs / d_ref, -1, 1)))
    for d_obs, d_ref in zip(d_obs_list, df_ds["d_ref_mm"].values)
]
labels = df_ds["label"].values

df_to_latex(
    pd.DataFrame({
        "狹縫代號": ["$\\varphi=\\arccos(d_{\\text{obs}}/d_{\\text{ref}})$ (\\unit{\\degree})"],
        **{str(lbl): [f"${phi:.2f}$"] for lbl, phi in zip(labels, phi_list)},
    }),
    caption="利用狹縫間距計算值反推狹縫與牆壁的（水平）夾角",
    path=REPORT,
)
