import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Parameters
AI, kI = 3.57, 0.964e-3
AE, kE = 3.35, 1.644e-3

# Intersection solver: A1(1 - e^{-k1P}) = A2(1 - e^{-k2P})
def f(P):
    return AI*(1 - np.exp(-kI*P)) - AE*(1 - np.exp(-kE*P))

# Solve for intersection (positive root)
P_int = fsolve(f, 1000)[0]
V_int = AI*(1 - np.exp(-kI*P_int))

# Generate data
P = np.linspace(0, P_int, 400)
VI = AI*(1 - np.exp(-kI*P))
VE = AE*(1 - np.exp(-kE*P))

# Plot
plt.plot(P, VI, label=f"Inhalation", color='blue')
plt.plot(P, VE, label=f"Exhalation", color='red')
# plt.plot(VI, P, label=f"Inhalation", color='blue')
# plt.plot(VE, P, label=f"Exhalation", color='red')

# Intersection point marker
plt.text(P_int, V_int*0.99, f"({P_int:.2f}, {V_int:.2f})", ha='center',va='top')
# plt.text(V_int*0.99, P_int, f"({V_int:.2f}, {P_int:.2f})", ha='right')

plt.xlabel("Pressure (Pa)")
plt.ylabel("Volume (L)")
# plt.ylabel("Pressure (Pa)")
# plt.xlabel("Volume (L)")

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title("V = A(1 - e$^{-kP})$")
plt.legend()
plt.show()