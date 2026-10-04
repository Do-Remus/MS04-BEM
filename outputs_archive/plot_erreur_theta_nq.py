import os

import matplotlib.pyplot as plt
import numpy as np

# 1. Charger les données
data_file = 'erreur_u_theta_vs_nq.txt'

save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

# Les valeurs de n_q sont écrites dans une ligne de commentaire du type :
# "# n_q_values: 1 2 3 4 5 6 8 10 15 20"
nq_values = []
with open(data_file, 'r') as f:
    for line in f:
        if line.startswith('# n_q_values:'):
            nq_values = [int(v) for v in line.split(':', 1)[1].split()]
            break

if not nq_values:
    raise RuntimeError("Impossible de trouver la ligne '# n_q_values:' dans " + data_file)

# np.loadtxt ignore automatiquement les lignes commençant par '#'
data = np.loadtxt(data_file)
theta = data[:, 0]
erreurs = data[:, 1:]  # une colonne par valeur de n_q

if erreurs.shape[1] != len(nq_values):
    raise RuntimeError("Nombre de colonnes d'erreur incohérent avec la liste des n_q")

# 2. Tracé : une courbe par n_q, superposées
fig, ax = plt.subplots(figsize=(9, 6))

cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(nq_values)))

for i, nq in enumerate(nq_values):
    ax.semilogy(theta, np.maximum(erreurs[:, i], 1e-18), color=colors[i],
                linewidth=1.5, label=f'$n_q$ = {nq}')

ax.set_xlabel(r'Angle $\theta$ (rad)', fontsize=13)
ax.set_ylabel(r'Erreur absolue $|u^+_{BEM} - u^+_{analytique}|$', fontsize=13)
ax.set_title(
    "Erreur de $u^+$ (BEM) vs solution analytique sur un cercle test\n"
    "en fonction de l'ordre de quadrature $n_q$ (pas de maillage fixe)",
    fontsize=13
)
ax.set_xlim(theta.min(), theta.max())
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9, ncol=2, loc='best', title=r'Ordre de quadrature')

fig.tight_layout()

png_path = os.path.join(save_dir, 'erreur_u_theta_vs_nq.png')
fig.savefig(png_path, dpi=200, bbox_inches='tight')
print('saved', png_path)

plt.show()