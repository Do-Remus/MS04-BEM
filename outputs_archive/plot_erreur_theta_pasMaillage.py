import os

import matplotlib.pyplot as plt
import numpy as np

# 1. Charger les données
data_file = 'erreur_u_theta_vs_pasMaillage.txt'

save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

# Les valeurs de pas de maillage sont écrites dans une ligne de commentaire du type :
# "# pasMaillage_values: 0.2 0.1 0.05 0.025 0.0125"
pas_values = []
with open(data_file, 'r') as f:
    for line in f:
        if line.startswith('# pasMaillage_values:'):
            pas_values = [float(v) for v in line.split(':', 1)[1].split()]
            break

if not pas_values:
    raise RuntimeError("Impossible de trouver la ligne '# pasMaillage_values:' dans " + data_file)

# np.loadtxt ignore automatiquement les lignes commençant par '#'
data = np.loadtxt(data_file)
theta = data[:, 0]
erreurs = data[:, 1:]  # une colonne par pas de maillage

if erreurs.shape[1] != len(pas_values):
    raise RuntimeError("Nombre de colonnes d'erreur incohérent avec la liste des pas de maillage")

# 2. Tracé : une courbe par pas de maillage, superposées
fig, ax = plt.subplots(figsize=(9, 6))

cmap = plt.cm.plasma
colors = cmap(np.linspace(0, 0.85, len(pas_values)))

for i, h in enumerate(pas_values):
    ax.semilogy(theta, np.maximum(erreurs[:, i], 1e-18), color=colors[i],
                linewidth=1.5, label=f'h = {h}')

ax.set_xlabel(r'Angle $\theta$ (rad)', fontsize=13)
ax.set_ylabel(r'Erreur absolue $|u^+_{BEM} - u^+_{analytique}|$', fontsize=13)
ax.set_title(
    "Erreur de $u^+$ (BEM) vs solution analytique sur un cercle test\n"
    "en fonction du pas de maillage de la frontière ($n_q$ fixe)",
    fontsize=13
)
ax.set_xlim(theta.min(), theta.max())
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9, loc='best', title='Pas de maillage')

fig.tight_layout()

png_path = os.path.join(save_dir, 'erreur_u_theta_vs_pasMaillage.png')
fig.savefig(png_path, dpi=200, bbox_inches='tight')
print('saved', png_path)

plt.show()