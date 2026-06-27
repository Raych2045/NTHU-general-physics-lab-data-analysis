import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import pm, typeA_uncertainty, df_to_latex

VERNIER_RES_CM = 0.005      # cm, vernier caliper least count (0.05 mm)
MICROMETER_RES_MM = 0.01    # mm, micrometer caliper least count
BALANCE_RES_G = 0.01        # g, electronic balance least count

BASE_DIR = Path(__file__).resolve().parent.parent
REPORT = BASE_DIR / "lab 1" / "report.tex"
REPORT.write_text("", encoding="utf-8")

DATA_PATH = BASE_DIR / "lab 1" / "data.xlsx"

df_pvc = pd.read_excel(DATA_PATH, sheet_name="pvc")
df_cylinder = pd.read_excel(DATA_PATH, sheet_name="cylinder")
df_steel_ball = pd.read_excel(DATA_PATH, sheet_name="steel_ball")


def stats(data, resolution):
    data = np.asarray(data, dtype=float)
    mean, uA = typeA_uncertainty(data)
    std = data.std(ddof=1)
    uB = resolution / np.sqrt(12)
    uC = np.sqrt(uA**2 + uB**2)
    return mean, std, uA, uB, uC


def stats_table(rows, caption):
    df_to_latex(
        pd.DataFrame({
            "Measurement item": [r[0] for r in rows],
            "Mean": [f"{r[1]:.5f}" for r in rows],
            "Std. dev.": [f"{r[2]:.3g}" for r in rows],
            "Type A unc.": [f"{r[3]:.3g}" for r in rows],
            "Type B unc.": [f"{r[4]:.3g}" for r in rows],
            "Combined unc.": [f"{r[5]:.3g}" for r in rows],
        }),
        caption=caption,
        path=REPORT,
    )


def qv_table(rows, caption):
    df_to_latex(
        pd.DataFrame({
            "Quantity": [r[0] for r in rows],
            "Value": [r[1] for r in rows],
        }),
        caption=caption,
        path=REPORT,
    )


def rel_propagate(value, terms):
    # terms: list of (power, x, ux); returns absolute uncertainty of value = f(...x...)
    return abs(value) * np.sqrt(sum((p * ux / x) ** 2 for p, x, ux in terms))


# =========================
# Table 1: PVC pipe -- raw measurement statistics
# =========================
inner_d = stats(df_pvc["inner_diameter_cm"], VERNIER_RES_CM)
outer_d = stats(df_pvc["outer_diameter_cm"], VERNIER_RES_CM)
length = stats(df_pvc["length_cm"], VERNIER_RES_CM)
mass_pvc = stats(df_pvc["mass_g"], BALANCE_RES_G)
buoyancy = stats(df_pvc["buoyancy_gw"], BALANCE_RES_G)

stats_table(
    [
        ("Inner diameter (\\unit{\\cm})", *inner_d),
        ("Outer diameter (\\unit{\\cm})", *outer_d),
        ("Length (\\unit{\\cm})", *length),
        ("Mass (\\unit{\\g})", *mass_pvc),
        ("Buoyancy reading (\\unit{gw})", *buoyancy),
    ],
    "PVC pipe -- summary statistics of the raw measurements",
)

# =========================
# Table 2: Cylinder -- raw measurement statistics
# =========================
diameter_cyl = stats(df_cylinder["diameter_cm"], VERNIER_RES_CM)
height_cyl = stats(df_cylinder["height_cm"], VERNIER_RES_CM)
mass_cyl = stats(df_cylinder["mass_g"], BALANCE_RES_G)

stats_table(
    [
        ("Diameter (\\unit{\\cm})", *diameter_cyl),
        ("Height (\\unit{\\cm})", *height_cyl),
        ("Mass (\\unit{\\g})", *mass_cyl),
    ],
    "Cylinder -- summary statistics of the raw measurements",
)

# =========================
# Table 3: Steel ball -- raw measurement statistics
# =========================
corrected_diameter_mm = (
    df_steel_ball["diameter_raw_mm"] - df_steel_ball["zero_reading_mm"]
)
diameter_ball = stats(corrected_diameter_mm, MICROMETER_RES_MM)
mass_ball = stats(df_steel_ball["mass_g"], BALANCE_RES_G)

stats_table(
    [
        ("Diameter (\\unit{\\mm})", *diameter_ball),
        ("Mass (\\unit{\\g})", *mass_ball),
    ],
    "Steel ball -- summary statistics of the raw measurements",
)

# =========================
# Tables 4+5 (merged): PVC pipe -- density calculation
# =========================
r_inner, u_r_inner = inner_d[0] / 2, inner_d[4] / 2
R_outer, u_R_outer = outer_d[0] / 2, outer_d[4] / 2

inner_area = np.pi * r_inner**2
u_inner_area = rel_propagate(inner_area, [(2, r_inner, u_r_inner)])

outer_area = np.pi * R_outer**2
u_outer_area = rel_propagate(outer_area, [(2, R_outer, u_R_outer)])

base_area_pvc = outer_area - inner_area
u_base_area_pvc = np.sqrt(u_outer_area**2 + u_inner_area**2)

volume_pvc = base_area_pvc * length[0]
u_volume_pvc = rel_propagate(
    volume_pvc, [(1, base_area_pvc, u_base_area_pvc), (1, length[0], length[4])]
)

density_pvc = mass_pvc[0] / volume_pvc
u_density_pvc = rel_propagate(
    density_pvc, [(1, mass_pvc[0], mass_pvc[4]), (1, volume_pvc, u_volume_pvc)]
)

density_pvc_arch = mass_pvc[0] / buoyancy[0]
u_density_pvc_arch = rel_propagate(
    density_pvc_arch, [(1, mass_pvc[0], mass_pvc[4]), (1, buoyancy[0], buoyancy[4])]
)

qv_table(
    [
        ("Inner radius (\\unit{\\cm})", pm(r_inner, u_r_inner)),
        ("Outer radius (\\unit{\\cm})", pm(R_outer, u_R_outer)),
        ("Inner circle area (\\unit{\\square\\cm})", pm(inner_area, u_inner_area)),
        ("Outer circle area (\\unit{\\square\\cm})", pm(outer_area, u_outer_area)),
        ("Base area (\\unit{\\square\\cm})", pm(base_area_pvc, u_base_area_pvc)),
        ("Volume (\\unit{\\cubic\\cm})", pm(volume_pvc, u_volume_pvc)),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_pvc, u_density_pvc)),
        ("Density via Archimedes' principle (\\unit{\\g/\\cubic\\cm})", pm(density_pvc_arch, u_density_pvc_arch)),
    ],
    "PVC pipe -- density calculation",
)

# =========================
# Table 6: Cylinder -- density calculation
# =========================
r_cyl, u_r_cyl = diameter_cyl[0] / 2, diameter_cyl[4] / 2

base_area_cyl = np.pi * r_cyl**2
u_base_area_cyl = rel_propagate(base_area_cyl, [(2, r_cyl, u_r_cyl)])

volume_cyl = base_area_cyl * height_cyl[0]
u_volume_cyl = rel_propagate(
    volume_cyl, [(1, base_area_cyl, u_base_area_cyl), (1, height_cyl[0], height_cyl[4])]
)

density_cyl = mass_cyl[0] / volume_cyl
u_density_cyl = rel_propagate(
    density_cyl, [(1, mass_cyl[0], mass_cyl[4]), (1, volume_cyl, u_volume_cyl)]
)

qv_table(
    [
        ("Radius (\\unit{\\cm})", pm(r_cyl, u_r_cyl)),
        ("Base area (\\unit{\\square\\cm})", pm(base_area_cyl, u_base_area_cyl)),
        ("Volume (\\unit{\\cubic\\cm})", pm(volume_cyl, u_volume_cyl)),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_cyl, u_density_cyl)),
    ],
    "Cylinder -- density calculation",
)

# =========================
# Table 7: Steel ball -- density calculation
# =========================
r_ball_mm, u_r_ball_mm = diameter_ball[0] / 2, diameter_ball[4] / 2
r_ball_cm, u_r_ball_cm = r_ball_mm / 10, u_r_ball_mm / 10

volume_ball = (4 / 3) * np.pi * r_ball_cm**3
u_volume_ball = rel_propagate(volume_ball, [(3, r_ball_cm, u_r_ball_cm)])

density_ball = mass_ball[0] / volume_ball
u_density_ball = rel_propagate(
    density_ball, [(1, mass_ball[0], mass_ball[4]), (1, volume_ball, u_volume_ball)]
)

qv_table(
    [
        ("Radius (\\unit{\\mm})", pm(r_ball_mm, u_r_ball_mm)),
        ("Volume (\\unit{\\cubic\\cm})", pm(volume_ball, u_volume_ball)),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_ball, u_density_ball)),
    ],
    "Steel ball -- density calculation",
)

# =========================
# Results summary (結果彙整): one Quantity & Value table per object
# =========================
qv_table(
    [
        ("Inner diameter (\\unit{\\cm})", pm(inner_d[0], inner_d[4])),
        ("Outer diameter (\\unit{\\cm})", pm(outer_d[0], outer_d[4])),
        ("Length (\\unit{\\cm})", pm(length[0], length[4])),
        ("Mass (\\unit{\\g})", pm(mass_pvc[0], mass_pvc[4])),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_pvc, u_density_pvc)),
        ("Buoyancy reading (\\unit{gw})", pm(buoyancy[0], buoyancy[4])),
        ("Density via Archimedes' principle (\\unit{\\g/\\cubic\\cm})", pm(density_pvc_arch, u_density_pvc_arch)),
    ],
    "Results summary -- PVC pipe",
)

qv_table(
    [
        ("Diameter (\\unit{\\cm})", pm(diameter_cyl[0], diameter_cyl[4])),
        ("Height (\\unit{\\cm})", pm(height_cyl[0], height_cyl[4])),
        ("Mass (\\unit{\\g})", pm(mass_cyl[0], mass_cyl[4])),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_cyl, u_density_cyl)),
    ],
    "Results summary -- Cylinder",
)

qv_table(
    [
        ("Diameter (\\unit{\\mm})", pm(diameter_ball[0], diameter_ball[4])),
        ("Mass (\\unit{\\g})", pm(mass_ball[0], mass_ball[4])),
        ("Density (\\unit{\\g/\\cubic\\cm})", pm(density_ball, u_density_ball)),
    ],
    "Results summary -- Steel ball",
)

# =========================
# Appendix: raw data tables
# =========================
df_to_latex(
    pd.DataFrame({
        "Inner diameter (\\unit{\\cm})": df_pvc["inner_diameter_cm"].map(lambda x: f"${x:.3f}$"),
        "Outer diameter (\\unit{\\cm})": df_pvc["outer_diameter_cm"].map(lambda x: f"${x:.3f}$"),
        "Length (\\unit{\\cm})": df_pvc["length_cm"].map(lambda x: f"${x:.3f}$"),
        "Mass (\\unit{\\g})": df_pvc["mass_g"].map(lambda x: f"${x:.2f}$"),
        "Buoyancy reading (\\unit{gw})": df_pvc["buoyancy_gw"].map(lambda x: f"${x:.2f}$"),
    }),
    caption="Raw measurements -- PVC pipe",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Diameter (\\unit{\\cm})": df_cylinder["diameter_cm"].map(lambda x: f"${x:.3f}$"),
        "Height (\\unit{\\cm})": df_cylinder["height_cm"].map(lambda x: f"${x:.3f}$"),
        "Mass (\\unit{\\g})": df_cylinder["mass_g"].map(lambda x: f"${x:.2f}$"),
    }),
    caption="Raw measurements -- Cylinder",
    path=REPORT,
)

df_to_latex(
    pd.DataFrame({
        "Zero reading (\\unit{\\mm})": df_steel_ball["zero_reading_mm"].map(lambda x: f"${x:.3f}$"),
        "Diameter reading (\\unit{\\mm})": df_steel_ball["diameter_raw_mm"].map(lambda x: f"${x:.3f}$"),
        "Corrected diameter (\\unit{\\mm})": corrected_diameter_mm.map(lambda x: f"${x:.3f}$"),
        "Mass (\\unit{\\g})": df_steel_ball["mass_g"].map(lambda x: f"${x:.2f}$"),
    }),
    caption="Raw measurements -- Steel ball",
    path=REPORT,
)
