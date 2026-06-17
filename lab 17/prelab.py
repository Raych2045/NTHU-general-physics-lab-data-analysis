import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

e_charge = 1.602176634e-19
c_light = 2.99792458e8
h_true = 6.62607015e-34
# lambda_blue = 460e-9
lambda_red = 650e-9

# =========================
# Helper functions
# =========================
def linear_func(x, a, b):
    return a * x + b

def compute_r_squared(x,a,b,y):
    y_predicted = linear_func(x, a, b)
    residuals = y - y_predicted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def weighted_linear_fit(x, y):
    popt, pcov = curve_fit(linear_func, x, y)
    slope, intercept = popt
    slope_err, intercept_err = np.sqrt(np.diag(pcov))
    cov_sb = pcov[0, 1]
    r_squared = compute_r_squared(x,slope,intercept,y)
    return slope, intercept, slope_err, intercept_err, r_squared, cov_sb

def find_tangent_region(V, I, window=5):
    dIdV = np.gradient(I, V)
    idx = np.argmax(dIdV)
    start = max(0, idx - window)
    end = min(len(V), idx + window)
    return start, end

def compute_relative_error(measured, true):
    return (measured - true) / true * 100.0

def pm(value, uncertainty):
    if uncertainty <= 0 or np.isnan(uncertainty):
        return f"${value}$"

    # ---------- helper: lab-style rounding ----------
    def round_uncertainty(u):
        exp = int(math.floor(math.log10(abs(u))))
        scale = 10 ** exp
        mant = u / scale  # between [1, 10)

        # keep 2 significant digits
        factor = 10
        mant_scaled = mant * factor

        next_digit = int((mant_scaled * 10) % 10)

        if next_digit == 0:
            mant_rounded = math.floor(mant_scaled) / factor
        else:
            mant_rounded = math.ceil(mant_scaled) / factor

        return mant_rounded * scale, exp

    # ---------- round uncertainty ----------
    u_rounded, u_exp = round_uncertainty(uncertainty)

    # decimal places for fixed-point rounding
    decimals = -u_exp + 1
    v_rounded = round(value, decimals)

    # ---------- decide scientific notation ----------
    use_sci = (decimals < 0) or (decimals > 3)

    if use_sci:
        sci_exp = int(math.floor(math.log10(abs(value))))
        v_scaled = v_rounded / (10 ** sci_exp)
        u_scaled = u_rounded / (10 ** sci_exp)

        fmt = f"{{:.{decimals+sci_exp}f}}" # |sci_exp|<=|decimals|
        return (
            f"$({fmt.format(v_scaled)} \\pm {fmt.format(u_scaled)})"
            f"\\times 10^{{{sci_exp}}}$"
        )

    else:
        fmt = f"{{:.{decimals}f}}"
        return f"${fmt.format(v_rounded)} \\pm {fmt.format(u_rounded)}$"

# for blue
# V = np.array([0.01, 1.01, 2.01, 2.51, 2.60, 2.65, 2.69, 2.73, 2.76, 2.79, 2.81, 2.84, 2.86])
# I = np.array([0, 0, 0, 0.242, 1.72, 3.20, 4.69, 6.18, 7.67, 9.17, 10.67, 12.17, 13.66])

# for red
V = np.array([0.012, 1.008, 1.697, 1.755, 1.810, 1.839, 1.860, 1.877, 1.892, 1.906, 1.919, 1.931, 1.945])
I = np.array([0, 0, 0.152, 0.606, 2.11, 3.60, 5.11, 6.61, 8.11, 9.61, 11.12, 12.62, 14.12])

idx_sort = np.argsort(V)
V = V[idx_sort]
I = I[idx_sort]
I=I*0.01

# Plot raw I-V
plt.figure(figsize=(9, 5))
plt.axhline(0,color="gray",alpha=0.7,linewidth=0.8)
plt.plot(V, I, "o-", label="Data points")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.title("I-V Curve of Red LED")

# Tangent method
start, end = find_tangent_region(V, I)
V_fit = V[start:end] # exclusive
I_fit = I[start:end]
print(V_fit)

s_tan, b_tan, s_tan_err, b_tan_err, r_tan, cov_sb = weighted_linear_fit(V_fit, I_fit)

# Turn-on voltage (I=0 => V = -b/s)
VF = -b_tan / s_tan

# Uncertainty propagation
VF_err = np.sqrt(
    (b_tan_err / s_tan)**2 +
    (b_tan * s_tan_err / s_tan**2)**2 -
    2 * (b_tan / s_tan**3) * cov_sb
)

# Plot tangent
V_line = np.linspace(VF, max(V), 100)
if b_tan<0:
    plt.plot(V_line, linear_func(V_line, s_tan, b_tan), color='black', linestyle="--", label=f'Tangent Fit Line: y = {s_tan:.3f} x - {abs(b_tan):.3f} (R squared = {r_tan:.4f})')
else:
    plt.plot(V_line, linear_func(V_line, s_tan, b_tan), color='black', linestyle="--", label=f'Tangent Fit Line: y = {s_tan:.3f} x + {b_tan:.3f} (R squared = {r_tan:.4f})')
plt.legend()
plt.show()

print(VF,VF_err,pm(VF,VF_err))

# =========================
# Planck constant
# =========================
K = e_charge * lambda_red / c_light
h = VF * K
h_err = VF_err * K

rel_err_h = compute_relative_error(h, h_true)
print(h,h_err,pm(h,h_err),rel_err_h)