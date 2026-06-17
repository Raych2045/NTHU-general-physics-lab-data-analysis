# Required libraries
import numpy as np
import matplotlib.pyplot as plt

# Constants
m = [1,2,3]
d = 10000 / 6
lambdas = [434, 410.2]  # wavelengths

# Theta range
theta_deg = np.linspace(-5, 5, 1000)
theta = np.deg2rad(theta_deg)

# Define lambda_m function
def lambda_m(m, theta, lam):
    argument = np.sin(theta) - m * lam / d
    argument = np.clip(argument, -1, 1)
    return -d/m * np.sin(np.arcsin(argument) - theta)

# Plot
plt.figure()

for M in m:
    for lam in lambdas:
        lambda_pos = lambda_m(M, theta, lam)
        lambda_neg = lambda_m(-M, theta, lam)
        Error = ((lambda_pos + lambda_neg)/(2*lam) - 1) * 100
        plt.plot(theta_deg, Error, label=f"m = {M}, λ = {lam} nm")

plt.xlabel(r"$\theta_{i}$ (degree)")
plt.ylabel("Error (%)")
plt.title(r"Error vs $\theta_{i}$ for Different Wavelengths")
plt.xlim(-5,5)
# plt.ylim(bottom=0)
plt.legend()
plt.show()