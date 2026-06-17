import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from package import df_to_latex

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "lab 21" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

REPORT = BASE_DIR / "lab 21" / "report.tex"
REPORT.write_text("", encoding="utf-8")

R = 0.105  # coil radius in meters


def B_single(x, center=0.0, direction=1):
    return direction * R**2 / (R**2 + (x - center)**2)**(3/2)


def simulate_field(subexp_id, x_m):
    if subexp_id == "A":
        return B_single(x_m)
    elif subexp_id == "B":
        d = 0.05
        return B_single(x_m, -d/2) + B_single(x_m, d/2)
    elif subexp_id == "C":
        d = 0.105
        return B_single(x_m, -d/2) + B_single(x_m, d/2)
    elif subexp_id == "D":
        d = 0.15
        return B_single(x_m, -d/2) + B_single(x_m, d/2)
    elif subexp_id == "E":
        d = 0.105
        return B_single(x_m, -d/2, direction=-1) + B_single(x_m, d/2)
    else:
        raise ValueError(f"Unknown subexp_id: {subexp_id}")


df = pd.read_excel(BASE_DIR / "lab 21" / "data.xlsx")

required_cols = {"subexp_id", "d", "Vpp"}
if not required_cols.issubset(df.columns):
    raise ValueError("Excel file must contain columns: subexp_id, d, Vpp")

for subexp_id in df["subexp_id"].unique():
    sub_df = df[df["subexp_id"] == subexp_id].copy()
    sub_df = sub_df.sort_values("d")

    x_cm = sub_df["d"].values
    Vpp = sub_df["Vpp"].values
    left, right = float(np.min(x_cm)), float(np.max(x_cm))

    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.plot(x_cm, Vpp, "o-", label="Data Point")
    plt.xlabel("x (cm)")
    plt.ylabel("Vpp (mV)")
    plt.xlim(left, right)
    plt.title(f"Subexperiment {subexp_id}")
    plt.legend()
    plt.savefig(FIGURES_DIR / f"data_subexp_{subexp_id}.png", dpi=300)
    plt.close()

    Vpp_norm = Vpp / np.max(np.abs(Vpp))
    x_sim = np.linspace(left / 100, right / 100, 1000)
    B_sim = simulate_field(subexp_id, x_sim)
    B_sim_norm = B_sim / np.max(np.abs(B_sim))

    plt.figure(figsize=(8, 5))
    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
    plt.scatter(x_cm, Vpp_norm, alpha=0.7, c="w", linewidths=1, edgecolors="r", label="Data Point (normalized)")
    plt.plot(x_sim * 100, B_sim_norm, label="Simulation (normalized)")
    plt.xlabel("x (cm)")
    plt.ylabel("Normalized amplitude")
    plt.xlim(left, right)
    if subexp_id == "E":
        plt.ylim(-1, 1)
    else:
        plt.ylim(0, 1)
    plt.title(f"Subexperiment {subexp_id}")
    plt.legend()
    plt.savefig(FIGURES_DIR / f"subexp_{subexp_id}.png", dpi=300)
    plt.close()

    table = pd.DataFrame({
        "$d$ (\\unit{\\cm})": sub_df["d"].map(lambda x: f"${x:.4g}$"),
        "Vpp (\\unit{\\mV})": sub_df["Vpp"].map(lambda x: f"${x:.4g}$"),
    })
    df_to_latex(table, path=REPORT, caption=f"Subexp {subexp_id}")

plt.figure(figsize=(8, 5))
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2, 3), useMathText=True)
for subexp_id in df["subexp_id"].unique():
    sub_df = df[df["subexp_id"] == subexp_id].copy()
    sub_df = sub_df.sort_values("d")
    x_cm = sub_df["d"].values
    Vpp = sub_df["Vpp"].values
    plt.plot(x_cm, np.abs(Vpp), ".-", label=f"Subexp {subexp_id}")
plt.xlabel("x (cm)")
plt.ylabel("Vpp (mV)")
plt.title("Data of all subexperiments")
plt.ylim(bottom=0)
plt.legend()
plt.savefig(FIGURES_DIR / "compare.png", dpi=300)
plt.close()
