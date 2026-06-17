import numpy as np
import matplotlib.pyplot as plt

# Given constants
m = 2
d = 10000 / 6 # nm
lam = 600  # wavelength (nm)

# Theta range in degrees
theta_deg = np.linspace(0, 10, 1000)
theta = np.deg2rad(theta_deg)  # convert to radians

# Define lambda_m function
def lambda_m(m, theta):
    argument = np.sin(theta) - m * lam / d
    
    # Avoid invalid arcsin values due to floating point rounding
    argument = np.clip(argument, -1, 1)
    
    return -d/m * np.sin(np.arcsin(argument) - theta)

# Compute lambda_m and lambda_-m
lambda_pos = lambda_m(m, theta)
lambda_neg = lambda_m(-m, theta)

# Compute Error
Error_m = ((lambda_pos + lambda_neg)/(2*lam) - 1) * 100

# Define lambda_m function
def theta_avg(m, theta):
    argument_pos = np.sin(theta) - m * lam / d
    # Avoid invalid arcsin values due to floating point rounding
    argument_pos = np.clip(argument_pos, -1, 1)

    argument_neg = np.sin(theta) + m * lam / d
    argument_neg = np.clip(argument_neg, -1, 1)
    argument_avg = (np.arcsin(argument_neg)-np.arcsin(argument_pos))/2
    return d/m * np.sin(argument_avg)

Error_avg = (theta_avg(m, theta)/lam - 1) * 100

# Plot
plt.figure()
plt.plot(theta_deg, Error_m, label=r"Averaging over $\lambda$")
plt.plot(theta_deg, Error_avg, label=r"Averaging over $\theta$")
plt.xlabel(r"$\theta_{i}$ (degree)")
plt.ylabel("Error (%)")
plt.xlim(0,10)
plt.ylim(bottom=0)
plt.legend()
plt.title(r"Error vs $\theta_{i}$")
plt.show()