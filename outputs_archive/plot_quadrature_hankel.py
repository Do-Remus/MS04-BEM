import os

import matplotlib.pyplot as plt
import numpy as np

# 1. Charger les données
data_file = 'convergence_quadrature_hankel.txt'
data = np.loadtxt(data_file)

save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

n = data[:, 0]
erreurs = data[:, 1:]  # une colonne par segment

# Labels des segments, dans l'ordre où le C++ les écrit
segments = [
    ('[0.01, 0.02]', 'champ proche'),
    ('[1.01, 1.02]', 'champ intermédiaire'),
    ('[100.01, 100.02]', 'champ lointain'),
]

colors = ['tab:red', 'tab:orange', 'tab:green']

# On plafonne les erreurs nulles pour l'échelle log
erreurs_plot = np.maximum(erreurs, 1e-18)

# 2. Tracé
fig, ax = plt.subplots(figsize=(8, 6))

for i, (label, distance) in enumerate(segments):
    if i < erreurs_plot.shape[1]:
        ax.semilogy(n, erreurs_plot[:, i], 'o-', color=colors[i], markersize=5,
                    label=f'{label}  ({distance})')

ax.set_xlabel('Ordre de quadrature n', fontsize=13)
ax.set_ylabel('Erreur absolue vs. référence (ordre élevé)', fontsize=13)
ax.set_title(
    "Convergence de la quadrature de Gauss-Legendre\n"
    "vers l'intégrale de $H_0^{(1)}(x)$ selon la distance à l'origine",
    fontsize=13
)
ax.set_xticks(np.arange(0, n.max() + 1, 2))
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=10, loc='best')

fig.tight_layout()

png_path = os.path.join(save_dir, 'convergence_quadrature_hankel.png')
fig.savefig(png_path, dpi=200, bbox_inches='tight')
print('saved', png_path)

plt.show()