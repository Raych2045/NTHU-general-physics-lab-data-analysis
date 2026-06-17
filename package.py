import numpy as np
import math

def linear_model(x, a, b):
    return a * x + b

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
    v_rounded = np.round(value, decimals)

    # ---------- decide scientific notation ----------
    sci_exp = int(math.floor(math.log10(abs(value))))
    use_sci = (sci_exp < -2) or (sci_exp > 2) or (decimals < 0)

    if use_sci:
        v_scaled = v_rounded / (10 ** sci_exp)
        u_scaled = u_rounded / (10 ** sci_exp)

        fmt = f"{{:.{decimals+sci_exp}f}}" # assume that |sci_exp|<=|decimals|
        return (
            f"$({fmt.format(v_scaled)} \\pm {fmt.format(u_scaled)})"
            f"\\times 10^{{{sci_exp}}}$"
        )

    else:
        fmt = f"{{:.{decimals}f}}"
        return f"${fmt.format(v_rounded)} \\pm {fmt.format(u_rounded)}$"

def percent(x, fmt=".2f"):
    return f"${x:{fmt}}\\%$"

def df_to_latex(df, caption, path):
    latex = df.to_latex(
        index=False,
        escape=False,
        column_format="c" * len(df.columns)
    )

    latex = latex.replace(
        "\\begin{tabular}",
        "\\begin{table}[htb!]\n\\centering\n"
        f"\\caption{{{caption}}}\n"
        "\\begin{tabular}"
    ).replace(
        "\\end{tabular}",
        "\\end{tabular}\n"
        "\\end{table}"
    )

    with path.open("a", encoding="utf-8") as f:
        f.write(latex)

def typeA_uncertainty(data):
    data = np.array(data, dtype=float)
    N = len(data)
    mean = np.mean(data)

    if N > 1:
        u = np.std(data, ddof=1) / np.sqrt(N)
    else:
        u = 0.0

    return mean, u

def compute_r_squared(x, a, b, y):
    y_predicted = linear_model(x, a, b)
    residuals = y - y_predicted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared