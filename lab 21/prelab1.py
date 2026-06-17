# On-axis magnetic field distribution for various coil configurations
import numpy as np
import matplotlib.pyplot as plt

# Physical constants and parameters
mu0 = 4 * np.pi * 1e-7
I = 2.0  # current in A
R = 0.105  # radius in meters (10.5 cm)
N = 200  # number of turns

# x range from -15 cm to 15 cm
x = np.linspace(-0.15, 0.15, 1000)

def B_single(x, center=0.0, direction=1):
    return direction * (mu0 * N * I * R**2) / (2 * (R**2 + (x - center)**2)**(3/2))

# (i) Single coil
B1 = B_single(x)

plt.figure(figsize=(8,5))
plt.plot(x * 100, B1)
plt.xlabel("x (cm)")
plt.ylabel("B (T)")
plt.xlim(-15,15)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
plt.title("On-axis Magnetic Field: Single Coil")
plt.show()

# (ii) Coil pair, same current direction
distances = [0.05, 0.105, 0.15]  # 5 cm, 10.5 cm, 15 cm in meters

plt.figure(figsize=(8,5))
for d in distances:
    B_pair = B_single(x, center=-d/2) + B_single(x, center=d/2)
    plt.plot(x * 100, B_pair, label=f"d = {d*100:.1f} cm")
plt.xlabel("x (cm)")
plt.ylabel("B (T)")
plt.xlim(-15,15)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
plt.legend()
plt.title("On-axis Magnetic Field: Coil Pair (Same Direction)")
plt.show()

# (iii) Coil pair, opposite current direction (d = 10.5 cm)
d = 0.105
B_opposite = B_single(x, center=d/2) + B_single(x, center=-d/2, direction=-1)

plt.figure(figsize=(8,5))
plt.plot(x * 100, B_opposite)
plt.xlabel("x (cm)")
plt.ylabel("B (T)")
plt.xlim(-15,15)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
plt.title("On-axis Magnetic Field: Coil Pair (Opposite Direction), d = 10.5 cm")
plt.show()