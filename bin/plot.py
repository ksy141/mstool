import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 10.0
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots(figsize=(4,3))
ax.plot()
ax.set_xlim([0,1])
ax.set_xlabel('xlabel')
ax.set_xticks([0.0000, 0.25000, 0.5000, 1.0000])
ax.set_xticklabels([0.0, 0.25, 0.5, 1.0])

ax.axhline()
ax.axvline()

