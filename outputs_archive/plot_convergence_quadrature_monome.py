import os

import matplotlib.pyplot as plt
import numpy as np

# 1. Charger les données
data_file = 'convergence_quadrature_monome.txt'
data = np.loadtxt(data_file)

save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

p = data[:, 0]
erreur = data[:, 1]

# Ordre de quadrature utilisé côté C++ (n=10) -> degré d'exactitude théorique = 2n-1
ordre_quadrature = 10
degre_exactitude = 2 * ordre_quadrature - 1

# On plafonne les erreurs nulles (ou proches de zéro machine) pour l'échelle log
erreur_plot = np.maximum(erreur, 1e-18)

# 2. Tracé
fig, ax = plt.subplots(figsize=(8, 6))

ax.semilogy(p, erreur_plot, 'o-', color='tab:blue', markersize=6,
            label=f'Erreur quadrature (n={ordre_quadrature})')

ax.axvline(degre_exactitude, color='tab:red', linestyle='--',
           label=f'Degré d\'exactitude théorique = {degre_exactitude}')

ax.axhline(1e-15, color='gray', linestyle=':', linewidth=1,
           label='Précision machine (~1e-15)')

ax.set_xlabel('Degré p du monôme $x^p$', fontsize=13)
ax.set_ylabel('Erreur absolue |quadrature - exacte|', fontsize=13)
ax.set_title(
    f'Erreur de la quadrature de Gauss-Legendre (n={ordre_quadrature}) sur [-1,1]\n'
    'pour les monômes de la base canonique',
    fontsize=13
)
ax.set_xticks(np.arange(0, p.max() + 1, 1))
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=10, loc='upper left')

fig.tight_layout()

png_path = os.path.join(save_dir, 'convergence_quadrature_monome.png')
fig.savefig(png_path, dpi=200, bbox_inches='tight')
print('saved', png_path)

plt.show()