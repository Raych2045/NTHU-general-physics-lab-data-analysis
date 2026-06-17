import numpy as np
import matplotlib.pyplot as plt

half_T = 0.5
total = 20 # 個半周期

# =======================
# Define piecewise function
# =======================
def y_piecewise(x_vals):
    y_vals = np.zeros_like(x_vals, dtype=float)

    # initial value
    y_prev = 0.0

    # process interval by interval
    for k in range(total):  # intervals [k, k+1], up to 10
        mask = (x_vals >= k*half_T) & (x_vals <= (k + 1)*half_T)

        x_segment = x_vals[mask]

        if k % 2 == 0:
            # even interval: growth
            y_segment = (1 - np.exp(-(x_segment - k*half_T))) + y_prev * np.exp(-(x_segment - k*half_T))
        else:
            # odd interval: decay
            y_segment = y_prev * np.exp(-(x_segment - k*half_T))

        y_vals[mask] = y_segment

        # update y_prev = y(k+1)
        if k % 2 == 0:
            y_prev = (1 - np.exp(-half_T)) + y_prev * np.exp(-half_T)
        else:
            y_prev = y_prev * np.exp(-half_T)

    return y_vals


# =======================
# Generate data
# =======================
x = np.linspace(0, total*half_T, 2000)
y = y_piecewise(x)

# =======================
# Plot
# =======================
plt.figure(figsize=(10,6))
plt.plot(x, y)
plt.xlabel("t/RC")
plt.ylabel(r"V/$\varepsilon_{0}$")
plt.xlim(0,total*half_T)
plt.ylim(0,1)
plt.title(f"Potential drop across the resistor versus time, using square wave with its period being {2*half_T}RC")
plt.show()