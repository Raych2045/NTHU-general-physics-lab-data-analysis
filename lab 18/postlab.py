from package import pm, percent
import numpy as np

# =========================
# Input values (EDIT HERE)
# =========================

w_central = 10.5   # central bright width
a = 0.08            # slit width
L = 778.5         # screen distance

# uncertainties
sigma_a = 0.01/np.sqrt(12)
sigma_L = 2.0/np.sqrt(12)

# =========================
# Derived uncertainty
# =========================

# =========================
# Wavelength
# =========================

lam = (w_central * a) / (2 * L)

# =========================
# Uncertainty propagation
# =========================

sigma_lam = lam * np.sqrt(
    (np.sqrt(0.01/6))**2 +
    (sigma_a / a)**2 +
    (sigma_L / L)**2
)

# =========================
# Output
# =========================

print("lambda =", lam)
print("sigma_lambda =", sigma_lam)
print(pm(lam,sigma_lam))
print("relative uncertainty =", percent(100*sigma_lam / lam))