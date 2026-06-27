"""
Lab 9 - Simple Harmonic Oscillation.

Four parts, 23 automated tables (-> report.tex) plus figures (-> figures/):

  Part 1  Spring constants
          static  (Hooke's law)  : Tables 1-3 (small), 4-6 (medium)
          dynamic (T^2 vs m)     : Tables 7-8 (small), 9-10 (medium)
          summary                : Table 11
  Part 2  Energy conservation of the vertical SHM   (figures only)
  Part 3  SHM under constant friction (cart turning points -> period and friction force)
          rail level             : Tables 12-14
          rail tilted            : Tables 15-16
  Part 4  Underdamped oscillation
          magnet high            : Table 17
          magnet low             : Table 18
          equal-mass blue cart   : Tables 19-20
          summary                : Table 21
  Discussion  relative-error comparisons             : Tables 22-23

All numbers are read from data.xlsx; only fixed apparatus constants are hard-coded.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import (pm, percent, df_to_latex, typeA_uncertainty,
                     linear_model, compute_r_squared)

BASE_DIR    = Path(__file__).resolve().parent.parent
LAB_DIR     = BASE_DIR / "lab 9"
FIGURES_DIR = LAB_DIR / "figures"
DATA        = LAB_DIR / "data.xlsx"
REPORT      = LAB_DIR / "report.tex"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT.write_text("", encoding="utf-8")          # clear on every run

g = 9.8                                            # gravitational acceleration

# ===========================================================================
# small helpers
# ===========================================================================
def kv_table(rows, caption):
    """Write a two-column (quantity / value) LaTeX table from a list of rows."""
    df = pd.DataFrame(rows, columns=["Quantity", "Value"])
    df_to_latex(df, caption=caption, path=REPORT)


def rel_error(obs, ref):
    """Signed relative error in percent."""
    return (obs - ref) / ref * 100.0


def setup_axes():
    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3),
                         useMathText=True)


def best_fit_label(s, b, r2):
    if b >= 0:
        return f'Best-fit line: y = {s:.4g} x + {b:.4g} (R² = {r2:.4f})'
    return f'Best-fit line: y = {s:.4g} x - {abs(b):.4g} (R² = {r2:.4f})'


# ===========================================================================
# PART 1 - spring constants
# ===========================================================================
def static_spring(sheet, k_sym, fig_name, title):
    """Static (Hooke's law) method -> returns (k, u_k)."""
    df = pd.read_excel(DATA, sheet_name=sheet)
    m_g, x_cm = df["m_g"].to_numpy(), df["x_cm"].to_numpy()

    # Table: raw hanging mass vs height above the floor
    raw = pd.DataFrame({
        r"Hanging mass $m$ (\unit{\g})":               [f"{v:.2f}" for v in m_g],
        r"Lower-coil height $x$ (\unit{\cm})":         [f"{v:.2f}" for v in x_cm],
    })
    df_to_latex(raw, caption="Height of the spring's lower coil above the "
                "floor for different hanging masses", path=REPORT)

    # extra restoring force and extension, relative to the first data point
    dF = (m_g[1:] - m_g[0]) / 1000.0 * g          # N
    dx = (x_cm[0] - x_cm[1:]) / 100.0             # m
    ext = pd.DataFrame({
        r"Restoring force $\Delta F=(m-m_{\text{ref}})g$ (\unit{\N})":
            [f"{v:.4f}" for v in dF],
        r"Extension $\Delta x=x_{\text{ref}}-x$ (\unit{\m})":
            [f"{v:.4f}" for v in dx],
    })
    df_to_latex(ext, caption="Additional restoring force and extension "
                "relative to the first measurement", path=REPORT)

    # linear fit  dF = k * dx
    (k, b), pcov = curve_fit(linear_model, dx, dF, p0=[10.0, 0.0])
    u_k, u_b = np.sqrt(np.diag(pcov))
    r2 = compute_r_squared(dx, k, b, dF)

    setup_axes()
    plt.scatter(dx, dF, label="Data Point")
    xf = np.linspace(dx.min(), dx.max(), 100)
    plt.plot(xf, k * xf + b, color='black', linestyle='--',
             label=best_fit_label(k, b, r2))
    plt.xlabel(r"$\Delta x$ (m)")
    plt.ylabel(r"$\Delta F$ (N)")
    plt.title(title)
    plt.legend()
    plt.savefig(FIGURES_DIR / fig_name, dpi=300)
    plt.close()

    # regression info  (R^2 omitted by request)
    kv_table([
        (rf"Slope ${k_sym}$ (\unit{{\N\per\m}})", pm(k, u_k)),
        (r"$y$-intercept (\unit{\N})",            pm(b, u_b)),
    ], caption=f"Linear-fit result for the static method")
    return k, u_k


def dynamic_spring(sheet, kd_sym, fig_name, title, ms_ref_g):
    """Dynamic (T^2 vs m) method -> returns (k_d, u_kd)."""
    df = pd.read_excel(DATA, sheet_name=sheet)
    m_g = df["m_g"].to_numpy()
    T1, T2 = df["T1_s"].to_numpy(), df["T2_s"].to_numpy()
    thirtyT = T2 - T1
    T = thirtyT / 30.0
    m_kg = m_g / 1000.0
    T2sq = T ** 2

    # Table: periods
    per = pd.DataFrame({
        r"Hanging mass $m$ (\unit{\g})":         [f"{v:.2f}" for v in m_g],
        r"$T_{1}$ (\unit{\s})":                  [f"{v:.2f}" for v in T1],
        r"$T_{2}$ (\unit{\s})":                  [f"{v:.2f}" for v in T2],
        r"$30T=T_{2}-T_{1}$ (\unit{\s})":        [f"{v:.2f}" for v in thirtyT],
        r"$T$ (\unit{\s})":                      [f"{v:.4g}" for v in T],
    })
    df_to_latex(per, caption="Oscillation period for different hanging "
                "masses ($T_{1}$: start of the 1st period, $T_{2}$: end of "
                "the 30th period)", path=REPORT)

    # linear fit  T^2 = s * m + b
    (s, b), pcov = curve_fit(linear_model, m_kg, T2sq, p0=[15.0, 0.0])
    u_s, u_b = np.sqrt(np.diag(pcov))
    r2 = compute_r_squared(m_kg, s, b, T2sq)

    kd = 4 * np.pi ** 2 / s
    u_kd = kd * u_s / s

    setup_axes()
    plt.scatter(m_kg, T2sq, label="Data Point")
    xf = np.linspace(m_kg.min(), m_kg.max(), 100)
    plt.plot(xf, s * xf + b, color='black', linestyle='--',
             label=best_fit_label(s, b, r2))
    plt.xlabel("m (kg)")
    plt.ylabel(r"$T^2$ (s$^2$)")
    plt.title(title)
    plt.legend()
    plt.savefig(FIGURES_DIR / fig_name, dpi=300)
    plt.close()

    # observed effective spring mass:  intercept b = s * m_s/3  ->  m_s = 3b/s
    # full propagation, including the slope-intercept covariance cov(s, b):
    #   var(m_s) = (3/s)^2 var_b + (3b/s^2)^2 var_s - 2(3/s)(3b/s^2) cov_sb
    ms_obs = 3 * b / s                                   # kg
    var_s, var_b, cov_sb = pcov[0, 0], pcov[1, 1], pcov[0, 1]
    u_ms = 3 * np.sqrt(var_b / s ** 2 + b ** 2 * var_s / s ** 4
                       - 2 * b * cov_sb / s ** 3)

    # regression info (R^2 omitted by request) merged with the m_s comparison;
    # the y-intercept b sits between the slope and k_d, since m_s,obs = 3b/s
    kv_table([
        (rf"Slope $s$ (\unit{{\s\squared\per\kg}})", pm(s, u_s)),
        (r"$y$-intercept $b$ (\unit{\s\squared})", pm(b, u_b)),
        (rf"${kd_sym}=4\pi^{{2}}/s$ (\unit{{\N\per\m}})", pm(kd, u_kd)),
        (r"$m_{s,\text{obs}}=3b/s$ (\unit{\g})", pm(ms_obs * 1000, u_ms * 1000)),
        (r"$m_{s,\text{ref}}$ (\unit{\g})",      f"${ms_ref_g:.2f}$"),
        ("Error", percent(rel_error(ms_obs * 1000, ms_ref_g), fmt="+.2f")),
    ], caption="Linear-fit result for the dynamic method and comparison of "
       "the observed and reference effective spring mass")
    return kd, u_kd


# ===========================================================================
# PART 2 - energy conservation (hard-coded constants, no error propagation)
# ===========================================================================
def energy_part(k_energy):
    df = pd.read_excel(DATA, sheet_name="energy_yt")
    t, y = df["t_s"].to_numpy(), df["y_m"].to_numpy()

    m   = 0.08036          # hanging mass (kg)
    m_s = 0.01316          # medium spring mass (kg)
    ydot = np.gradient(y, t)
    KE = 0.5 * (m + m_s / 3) * ydot ** 2
    PE = 0.5 * k_energy * y ** 2

    # y - t curve
    setup_axes()
    plt.plot(t, y, color="tab:blue")
    plt.xlabel("Time t (s)")
    plt.ylabel("Displacement y (m)")
    plt.xlim(left=0)
    plt.title("Vertical simple harmonic oscillation")
    plt.savefig(FIGURES_DIR / "y-t.png", dpi=300)
    plt.close()

    # kinetic / potential energy line plot
    setup_axes()
    plt.plot(t, KE, label="Kinetic energy")
    plt.plot(t, PE, label="Potential energy")
    plt.xlabel("Time t (s)")
    plt.ylabel("Energy E (J)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.title("Kinetic and potential energy versus time")
    plt.legend()
    plt.savefig(FIGURES_DIR / "energy-t.png", dpi=300)
    plt.close()

    # stacked area chart (total mechanical energy)
    setup_axes()
    plt.stackplot(t, KE, PE, labels=["Kinetic energy", "Potential energy"])
    plt.xlabel("Time t (s)")
    plt.ylabel("Energy E (J)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.title("Total mechanical energy versus time")
    plt.legend(loc="upper right")
    plt.savefig(FIGURES_DIR / "conservation.png", dpi=300)
    plt.close()


# ===========================================================================
# PART 3 / 4 - turning-point detection (cart oscillation under friction)
# ===========================================================================
def turning_points(t, x):
    """Return (t, x) at the alternating maxima/minima of the trace."""
    rng = x.max() - x.min()
    dt = float(np.median(np.diff(t)))              # Median of sampling step size
    distance = max(1, int(0.3 / dt))               # >= 0.3 s apart
    prominence = 0.1 * rng
    pk, _ = find_peaks(x,  prominence=prominence, distance=distance)
    tr, _ = find_peaks(-x, prominence=prominence, distance=distance)
    idx = np.sort(np.concatenate([pk, tr]))
    return t[idx], x[idx]


def friction_case(sheet, fig_name, fig_title, turn_caption):
    """Detect turning points, write the turning-point table, and return the
    period T and the per-cycle amplitude drop D (both with uncertainty)."""
    df = pd.read_excel(DATA, sheet_name=sheet)
    t, x = df["t_s"].to_numpy(), df["x_m"].to_numpy()
    tt, xt = turning_points(t, x)

    # x - t curve with the detected turning points marked
    setup_axes()
    plt.plot(t, x, color="tab:blue", label="$x(t)$")
    plt.scatter(tt, xt, marker="x", color="red", zorder=3, s=20, label="Turning point")
    plt.xlabel("Time t (s)")
    plt.ylabel("Displacement x (m)")
    plt.xlim(left=0)
    plt.title(fig_title)
    plt.legend()
    plt.savefig(FIGURES_DIR / fig_name, dpi=300)
    plt.close()

    # turning-point table
    tab = pd.DataFrame({
        r"$t$ (\unit{\s})": [f"{v:.4f}" for v in tt],
        r"$x$ (\unit{\m})": [f"{v:.5f}" for v in xt],
    })
    df_to_latex(tab, caption=turn_caption, path=REPORT)

    # period: each gap between successive turning points is T/2
    half = np.diff(tt)
    _, u_half = typeA_uncertainty(half)
    T, u_T = 2 * np.mean(half), 2 * u_half

    # amplitude drop over one full period: distance between turning points
    # spaced one entry apart (peak->peak, trough->trough)
    drops = np.abs(xt[2:] - xt[:-2])
    D, u_D = typeA_uncertainty(drops)
    return (T, u_T), (D, u_D)


def theoretical_period(M, m_eff, k_eff, u_keff):
    """T_theo = 2 pi sqrt((M + m_eff/3)/k_eff)  (masses in kg)."""
    eff = M + m_eff / 3
    T = 2 * np.pi * np.sqrt(eff / k_eff)
    u_T = np.pi * np.sqrt(eff / k_eff ** 3) * u_keff
    return T, u_T


# ===========================================================================
# PART 4 - underdamped fit
# ===========================================================================
def damped_model(t, A, B, C, D, E):
    return A * np.exp(E * t) * np.cos(B * t + C) + D


def damped_case(sheet, fig_name, fig_title, p0, k, u_k):
    """Fit x(t)=A e^{Et} cos(Bt+C)+D and propagate to omega0, zeta, M, c, T.
    Returns a dict of (value, uncertainty) pairs plus the raw fit parameters."""
    df = pd.read_excel(DATA, sheet_name=sheet)
    t, x = df["t_s"].to_numpy(), df["x_m"].to_numpy()
    popt, pcov = curve_fit(damped_model, t, x, p0=p0, maxfev=20000)
    A, B, C, D, E = popt
    perr = np.sqrt(np.diag(pcov))

    setup_axes()
    plt.scatter(t, x, s=8, label="Data Point")
    tf = np.linspace(t.min(), t.max(), 1000)
    plt.plot(tf, damped_model(tf, *popt), color='black',
             label="Best-fit curve")
    plt.xlabel("Time t (s)")
    plt.ylabel("Displacement x (m)")
    plt.xlim(left=0)
    plt.title(fig_title)
    plt.legend()
    plt.savefig(FIGURES_DIR / fig_name, dpi=300)
    plt.close()

    # covariance entries needed for B and E
    var_B, var_E, cov_BE = pcov[1, 1], pcov[4, 4], pcov[1, 4]

    # omega0^2 = B^2 + E^2  (kept for the M / c calculations below)
    omega = B ** 2 + E ** 2
    u_omega = np.sqrt(4 * B ** 2 * var_B + 4 * E ** 2 * var_E
                      + 8 * B * E * cov_BE)

    # natural angular frequency omega0 = sqrt(B^2 + E^2)
    omega0 = np.sqrt(omega)
    u_omega0 = u_omega / (2 * omega0)

    # damping ratio zeta = sqrt(E^2 / omega)
    zeta = np.sqrt(E ** 2 / omega)
    dz_dB = -E ** 2 * B / zeta / omega ** 2
    dz_dE = B ** 2 * E / zeta / omega ** 2
    u_zeta = np.sqrt(dz_dB ** 2 * var_B + dz_dE ** 2 * var_E
                     + 2 * dz_dB * dz_dE * cov_BE)

    # mass M = k / omega
    M = k / omega
    u_M = M * np.sqrt((u_k / k) ** 2 + (u_omega / omega) ** 2)

    # damping coefficient c = -2 M E = -2 k E / omega
    c = -2 * M * E
    dc_dk = -2 * E / omega
    dc_dB = 4 * k * B * E / omega ** 2
    dc_dE = 2 * k * (E ** 2 - B ** 2) / omega ** 2
    u_c = np.sqrt(dc_dB ** 2 * var_B + dc_dE ** 2 * var_E
                  + dc_dk ** 2 * u_k ** 2 + 2 * dc_dB * dc_dE * cov_BE)

    # quasi-period T = 2 pi / B
    T = 2 * np.pi / B
    u_T = 2 * np.pi / B ** 2 * np.sqrt(var_B)

    return {
        "A": (A, perr[0]), "B": (B, perr[1]), "C": (C, perr[2]),
        "D": (D, perr[3]), "E": (E, perr[4]),
        "omega": (omega, u_omega), "omega0": (omega0, u_omega0),
        "zeta": (zeta, u_zeta),
        "M": (M, u_M), "c": (c, u_c), "T": (T, u_T),
    }


# ===========================================================================
# main flow  (writes tables to report.tex in document order 1..23)
# ===========================================================================
# ---- Part 1: spring constants ---------------------------------------------
k_s1, u_ks1 = static_spring("static_small",  "k_{s,1}", "small ks.png",
                            "Hooke's law fit, small-k spring")
k_s2, u_ks2 = static_spring("static_medium", "k_{s,2}", "medium ks.png",
                            "Hooke's law fit, medium-k spring")
k_d1, u_kd1 = dynamic_spring("dynamic_small", "k_{d,1}",
                             "small kd.png", "Squared period vs mass, "
                             "small-k spring", ms_ref_g=10.40)
k_d2, u_kd2 = dynamic_spring("dynamic_medium", "k_{d,2}",
                             "medium kd.png", "Squared period vs mass, "
                             "medium-k spring", ms_ref_g=13.16)

# Table 13 - summary of the four spring-constant estimates
kv_table([
    (r"$k_{s,1}$ (\unit{\N\per\m})", pm(k_s1, u_ks1)),
    (r"$k_{d,1}$ (\unit{\N\per\m})", pm(k_d1, u_kd1)),
    (r"$k_{s,2}$ (\unit{\N\per\m})", pm(k_s2, u_ks2)),
    (r"$k_{d,2}$ (\unit{\N\per\m})", pm(k_d2, u_kd2)),
], caption="Comparison of static and dynamic spring constants")

# ---- choose which spring-constant estimates feed k_eff --------------------
# Two springs on both sides of the cart: one low-k, one medium-k.  Defaults reproduce
# the report (k_eff = k_d1 + k_s2).  Swap to k_s1 / k_d2 here if desired.
K_SMALL,  U_SMALL  = k_d1, u_kd1
K_MEDIUM, U_MEDIUM = k_s2, u_ks2
k_eff_two  = K_SMALL + K_MEDIUM
u_eff_two  = np.hypot(U_SMALL, U_MEDIUM)
k_eff_dbl  = 2 * k_s2                       # two medium springs (Part 4 / blue)
u_eff_dbl  = 2 * u_ks2

# ---- Part 2: energy conservation (uses k_s2 best estimate) ----------------
energy_part(k_energy=k_s2)

# ---- Part 3: constant friction --------------------------------------------
M_BLUE_LEVEL = 0.58718      # kg
M_BLUE_TILT  = 0.59986
M_EFF_TWO    = 0.02338      # two-spring total mass (kg)
THETA_DEG    = 43.7631      # rail tilt angle (delta theta neglected)

# acceleration of the unloaded cart (shared by Tables 16, 18, 22)
acc = pd.read_excel(DATA, sheet_name="accel")
a_all = np.concatenate([acc["right_cm_s2"].to_numpy(),
                        acc["left_cm_s2"].to_numpy()]) / 100.0   # m/s^2
a_mean, u_a = typeA_uncertainty(a_all)
abs_a = abs(a_mean)

# -- rail level (Tables 14, 15, 16)
(T_lv, uT_lv), (D_lv, uD_lv) = friction_case(
    "friction_level", "level.png",
    "SHM under constant friction (rail level)",
    "Turning times and turning points of the whole motion (rail level)")

# Table 15 - acceleration measurements
acc_tab = pd.DataFrame({
    r"Cart rightward (\unit{\cm\per\s\squared})":
        [f"{v:.2f}" for v in acc["right_cm_s2"]],
    r"Cart leftward (\unit{\cm\per\s\squared})":
        [f"{v:.2f}" for v in acc["left_cm_s2"]],
})
df_to_latex(acc_tab, caption="Measured acceleration of the unloaded cart",
            path=REPORT)

Tth_lv, uTth_lv = theoretical_period(M_BLUE_LEVEL, M_EFF_TWO, k_eff_two, u_eff_two)
f_lv  = k_eff_two * D_lv / 4
uf_lv = 0.25 * np.hypot(k_eff_two * uD_lv, D_lv * u_eff_two)
fp    = M_BLUE_LEVEL * abs_a
ufp   = M_BLUE_LEVEL * u_a
kv_table([
    (r"Period $T_{\text{obs}}$ (\unit{\s})", pm(T_lv, uT_lv)),
    (r"$M_{\text{blue}}$ (\unit{\g})", r"$587.18$"),
    (r"$m_{\text{eff}}$ (\unit{\g})", r"$23.38$"),
    (r"$k_{\text{eff}}=k_{d,1}+k_{s,2}$ (\unit{\N\per\m})", pm(k_eff_two, u_eff_two)),
    (r"$T_{\text{theo}}=2\pi \sqrt{ \dfrac{M+m_{\text{eff}} /3}{k_{\text{eff}}} }$ (\unit{\s})", pm(Tth_lv, uTth_lv)),
    (r"Error of $T_{\text{obs}}$", percent(rel_error(T_lv, Tth_lv), fmt="+.2f")),
    (r"$D$ (\unit{\m})", pm(D_lv, uD_lv)),
    (r"$f=k_{\text{eff}}D/4$ (\unit{\N})", pm(f_lv, uf_lv)),
    (r"$a$ (\unit{\m\per\s\squared})", pm(a_mean, u_a)),
    (r"$f'=M_{\text{blue}}|a|$ (\unit{\N})", pm(fp, ufp)),
    (r"Error of $f$ relative to $f'$", percent(rel_error(f_lv, fp), fmt="+.2f")),
], caption="Period and friction-force analysis with the rail level")

# -- rail tilted (Tables 17, 18)
(T_tl, uT_tl), (D_tl, uD_tl) = friction_case(
    "friction_tilt", "tilt.png",
    "SHM under a smaller constant friction force (rail tilted)",
    "Turning times and turning points of the whole motion (rail tilted)")

Tth_tl, uTth_tl = theoretical_period(M_BLUE_TILT, M_EFF_TWO, k_eff_two, u_eff_two)
f_tl  = k_eff_two * D_tl / 4
uf_tl = 0.25 * np.hypot(k_eff_two * uD_tl, D_tl * u_eff_two)
cos_th = np.cos(np.radians(THETA_DEG))
fpp    = M_BLUE_TILT * abs_a * cos_th
ufpp   = M_BLUE_TILT * cos_th * u_a
kv_table([
    (r"Period $T_{\text{obs}}$ (\unit{\s})", pm(T_tl, uT_tl)),
    (r"$M_{\text{blue}}'$ (\unit{\g})", r"$599.86$"),
    (r"$T_{\text{theo}}$ (\unit{\s})", pm(Tth_tl, uTth_tl)),
    (r"Error of $T_{\text{obs}}$", percent(rel_error(T_tl, Tth_tl), fmt="+.2f")),
    (r"$D$ (\unit{\m})", pm(D_tl, uD_tl)),
    (r"$f=k_{\text{eff}}D/4$ (\unit{\N})", pm(f_tl, uf_tl)),
    (r"$\theta$ (\unit{\degree})", rf"${THETA_DEG:.4f}$"),
    (r"$f''=M_{\text{blue}}'|a|\cos\theta$ (\unit{\N})", pm(fpp, ufpp)),
    (r"Error of $f$ relative to $f''$", percent(rel_error(f_tl, fpp), fmt="+.2f")),
], caption="Period and friction-force analysis with the rail tilted")

# ---- Part 4: underdamped oscillation --------------------------------------
M_RED   = 0.44664           # kg
M_EFF_D = 0.02632           # two medium springs (kg)
M_eff_total_g = (M_RED + M_EFF_D / 3) * 1000     # M_eff = M_red + m_eff/3

# Table 19 - magnet high
hi = damped_case("damped_high", "high.png", "Damped oscillation (magnet high)",
                 p0=[0.10049, 4.9031, -2.8264, 0.002297, -0.5969],
                 k=k_eff_dbl, u_k=u_eff_dbl)
M_hi_g = hi["M"][0] * 1000
kv_table([
    (r"$m_{\text{eff}}$ (\unit{\g})", r"$26.32$"),
    (r"$k_{\text{eff}}=2k_{s,2}$ (\unit{\N\per\m})", pm(k_eff_dbl, u_eff_dbl)),
    (r"$M_{\text{red}}$ (\unit{\g})", r"$446.64$"),
    (r"$M_{\text{eff}}=M_{\text{red}}+m_{\text{eff}}/3$ (\unit{\g})",
        rf"${M_eff_total_g:.2f}$"),
    (r"$A$", pm(*hi["A"])),
    (r"$B$", pm(*hi["B"])),
    (r"$C$", pm(*hi["C"])),
    (r"$D$", pm(*hi["D"])),
    (r"$E$", pm(*hi["E"])),
    (r"$\omega_{0}$ (\unit{\radian\per\s})", pm(*hi["omega0"])),
    (r"Damping ratio $\zeta$", pm(*hi["zeta"]) + r" ($<1$)"),
    (r"$M$ (\unit{\g})", pm(M_hi_g, hi["M"][1] * 1000)),
    (r"Error of $M$ relative to $M_{\text{eff}}$",
        percent(rel_error(M_hi_g, M_eff_total_g), fmt="+.2f")),
    (r"Damping coefficient $c$ (\unit{\kg\per\s})", pm(*hi["c"])),
    (r"Quasi-period $T$ (\unit{\s})", pm(*hi["T"])),
], caption="Curve-fitting and derived results, magnet at the higher position")

# Table 20 - magnet low
lo = damped_case("damped_low", "low.png", "Damped oscillation (magnet low)",
                 p0=[0.08526, 4.8264, -2.7911, -0.002234, -1.1118],
                 k=k_eff_dbl, u_k=u_eff_dbl)
M_lo_g = lo["M"][0] * 1000
kv_table([
    (r"$A$", pm(*lo["A"])),
    (r"$B$", pm(*lo["B"])),
    (r"$C$", pm(*lo["C"])),
    (r"$D$", pm(*lo["D"])),
    (r"$E$", pm(*lo["E"])),
    (r"$\omega_{0}$ (\unit{\radian\per\s})", pm(*lo["omega0"])),
    (r"Damping ratio $\zeta$", pm(*lo["zeta"]) + r" ($<1$)"),
    (r"$M$ (\unit{\g})", pm(M_lo_g, lo["M"][1] * 1000)),
    (r"Error of $M$ relative to $M_{\text{eff}}$",
        percent(rel_error(M_lo_g, M_eff_total_g), fmt="+.2f")),
    (r"Damping coefficient $c$ (\unit{\kg\per\s})", pm(*lo["c"])),
    (r"Quasi-period $T$ (\unit{\s})", pm(*lo["T"])),
], caption="Curve-fitting and derived results, magnet at the lower position")

# -- equal-mass blue cart (Tables 21, 22) -- uses the two-medium-spring k_eff
M_BLUE_EQ = 0.44776         # kg
(T_bl, uT_bl), (D_bl, uD_bl) = friction_case(
    "friction_blue", "blue cart.png",
    "SHM of the equal-mass blue cart under constant friction",
    "Turning times and turning points of the whole motion (equal-mass blue cart)")

Tth_bl, uTth_bl = theoretical_period(M_BLUE_EQ, M_EFF_D, k_eff_dbl, u_eff_dbl)
f_bl  = k_eff_dbl * D_bl / 4
uf_bl = 0.25 * np.hypot(k_eff_dbl * uD_bl, D_bl * u_eff_dbl)
fppp  = M_BLUE_EQ * abs_a
ufppp = M_BLUE_EQ * u_a
kv_table([
    (r"Period $T_{\text{obs}}$ (\unit{\s})", pm(T_bl, uT_bl)),
    (r"$M_{\text{blue}}''$ (\unit{\g})", r"$447.76$"),
    (r"$T_{\text{theo}}$ (\unit{\s})", pm(Tth_bl, uTth_bl)),
    (r"Error of $T_{\text{obs}}$", percent(rel_error(T_bl, Tth_bl), fmt="+.2f")),
    (r"$D$ (\unit{\m})", pm(D_bl, uD_bl)),
    (r"$f=k_{\text{eff}}D/4$ (\unit{\N})", pm(f_bl, uf_bl)),
    (r"$f'''=M_{\text{blue}}''|a|$ (\unit{\N})", pm(fppp, ufppp)),
    (r"Error of $f$ relative to $f'''$", percent(rel_error(f_bl, fppp), fmt="+.2f")),
], caption="Period and friction-force analysis for the equal-mass blue cart")

# Table 23 - (quasi-)period comparison
CONFIGS = ["Magnet at the lower position",
           "Magnet at the higher position",
           "No magnet (equal-mass blue cart)"]
period_df = pd.DataFrame({
    "Configuration": CONFIGS,
    r"Quasi-period (\unit{\s})": [pm(*lo["T"]), pm(*hi["T"]), pm(T_bl, uT_bl)],
})
df_to_latex(period_df,
            caption="Quasi-period of the cart in different configurations",
            path=REPORT)

# ---- discussion summary tables --------------------------------------------
FRICTION_CONFIGS = ["Rail level",
                    "Rail tilted",
                    "Equal-mass blue cart (no magnet)"]

# Table 24 - relative error of the friction force
fric_err_df = pd.DataFrame({
    "Configuration": FRICTION_CONFIGS,
    r"Relative error of $f$": [
        percent(rel_error(f_lv, fp),   fmt="+.2f"),
        percent(rel_error(f_tl, fpp),  fmt="+.2f"),
        percent(rel_error(f_bl, fppp), fmt="+.2f"),
    ],
})
df_to_latex(fric_err_df, caption="Comparison of the relative error of the "
            "friction force under different conditions", path=REPORT)

# Table 25 - relative error of the oscillation period
period_err_df = pd.DataFrame({
    "Configuration": FRICTION_CONFIGS,
    r"Relative error of $T_{\text{obs}}$": [
        percent(rel_error(T_lv, Tth_lv), fmt="+.2f"),
        percent(rel_error(T_tl, Tth_tl), fmt="+.2f"),
        percent(rel_error(T_bl, Tth_bl), fmt="+.2f"),
    ],
})
df_to_latex(period_err_df, caption="Comparison of the relative error of the "
            "oscillation period under different conditions", path=REPORT)
