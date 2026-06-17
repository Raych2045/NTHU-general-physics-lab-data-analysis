import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# ============================================================
# Parameters
# ============================================================

# Signal generator
V0 = 1.0          # signal amplitude [V]
Rg = 600.0         # generator internal resistance [ohm]

# Channel 1 resistor
R1 = 100.0        # [ohm]

# Driving (multi-turn) coil
Rd = 7.2         # resistance [ohm]
Ld = 16.3e-3         # inductance [H]

# Pick-up coil RLC circuit
Rp = 25.0         # resistance [ohm]
Lp = 10e-3         # inductance [H]
Cp = 0.1e-6         # capacitance [F]

# Mutual inductance
M = 1e-4          # [H]

# ============================================================
# Frequency range
# ============================================================

f = np.linspace(2000, 8000, 5000)   # frequency [Hz]
w = 2 * np.pi * f                   # angular frequency

# ============================================================
# Driving circuit impedance
# ============================================================

Zd = (Rg + R1 + Rd) + 1j * w * Ld

# ============================================================
# Pick-up circuit impedance
# ============================================================

Zp = Rp + 1j * (w * Lp - 1 / (w * Cp))

# ============================================================
# Channel 2 voltage
# ============================================================

Vch1 = -(M * V0 / Cp) / (Zd * Zp)

# Amplitude
Vch1_amp = np.abs(Vch1)

# ============================================================
# Numerical resonant frequency
# ============================================================

peak_index = np.argmax(Vch1_amp)
f_res_numeric = f[peak_index]

# Theoretical resonance
f_res_theory = np.sqrt(1/(Lp * Cp)- Rp**2 /(2 * Lp**2)) / (2 * np.pi)

print(f"Numerical resonant frequency : {f_res_numeric:.2f} Hz")
print(f"Theoretical resonant frequency: {f_res_theory:.2f} Hz")

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(f, Vch1_amp)

plt.axvline(
    f_res_numeric,
    linestyle='--',
    label=f"Numerical resonance = {f_res_numeric:.2f} Hz"
)

plt.xlabel("Frequency (Hz)")
plt.ylabel(r"|V$_{\text{CH1}}$| (V)")
plt.title(r"Frequency Response of |V$_{\text{CH1}}$|")
plt.grid(True)
plt.ylim(bottom=0)
plt.xlim(2000,8000)
plt.legend()
plt.tight_layout()
plt.show()