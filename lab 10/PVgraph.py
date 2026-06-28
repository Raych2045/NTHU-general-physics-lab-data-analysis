import matplotlib.pyplot as plt

# Data
# low
# T = [302.15,302.15,326.15,326.15,302.15]
# P = [100.09,103.99,103.99,99.96,99.96]
# V = [184.6,176.5,185.5,192,183]
# high
T = [302.15,302.15,358.15,358.15,302.15]
P = [100.21,103.99,103.99,100.21,100.21]
V = [182.5,176,196.9,202.9,182]

# Determine high and low temperature groups
low_T = min(T)
high_T = max(T)

colors = ["red" if t == high_T else "blue" for t in T]

# Scatter points
plt.scatter(V, P, c=colors)

# Coordinate labels 

plt.text(V[0], P[0]+0.05, f"A ({V[0]}, {P[0]})", ha='left')
plt.text(V[1], P[1]+0.05, f"B ({V[1]}, {P[1]})", ha='left')
plt.text(V[2], P[2]+0.05, f"C ({V[2]}, {P[2]})", ha='left')
plt.text(V[3], P[3]+0.05, f"D ({V[3]}, {P[3]})", ha='left')

plt.text(V[4], P[4]+0.05, f"A' ({V[4]}, {P[4]})", ha='right')

# Arrows using annotate
for i in range(len(V)-1):
    plt.annotate("",
                 xy=(V[i+1], P[i+1]),
                 xytext=(V[i], P[i]),
                 arrowprops=dict(arrowstyle='->',
                connectionstyle='arc3, rad=0',
                color='gray',
                shrinkA=7, shrinkB=7))

# Legend (red = high T, blue = low T)
plt.scatter([], [], color="red", label=f"Higher T ({high_T} K)")
plt.scatter([], [], color="blue", label=f"Lower T ({low_T} K)")
plt.legend()

plt.xlabel("Volume (mL)")
plt.ylabel("Pressure (kPa)")
plt.title("PV Diagram with Arrows and Temperature Coloring")

plt.tight_layout()
plt.show()