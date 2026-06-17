import matplotlib.pyplot as plt
import matplotlib as mpl

# Enable LaTeX for all text
mpl.rcParams['text.usetex'] = True

# Optional: Set default font to serif (classic Computer Modern font)
mpl.rcParams['font.family'] = 'serif'

# Optional: Add specific LaTeX packages to the preamble
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'

# Now you can create your plots, and all text will be rendered with LaTeX
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Signal ($S(t)$)')
plt.title(r'An equation: $\frac{1}{\sqrt{\pi t_H}} \exp\left(-\frac{t^2}{t_H^2}\right)$')
plt.show()