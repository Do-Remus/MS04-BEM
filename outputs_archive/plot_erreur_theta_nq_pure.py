import os

import matplotlib.pyplot as plt
import numpy as np

# 1. Charger les données
data_file = 'erreur_u_theta_vs_nq_quadrature_pure.txt'

save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

# Les valeurs de n_q sont écrites dans une ligne de commentaire du type :
# "# n_q_values: 1 2 3 4 5 6 8 10 15 20"
nq_values = []
nq_reference = None
with open(data_file, 'r') as f:
    for line in f:
        if line.startswith('# n_q_values:'):
            nq_values = [int(v) for v in line.split(':', 1)[1].split()]
        if 'nqReference=' in line:
            nq_reference = line.split('nqReference=', 1)[1].split()[0]
        if nq_values and nq_reference is not None:
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
ax.set_ylabel(r'Erreur absolue $|u^+_{n_q} - u^+_{ref}|$', fontsize=13)
ref_str = f' (réf. $n_q$={nq_reference})' if nq_reference else ''
ax.set_title(
    "Erreur de quadrature PURE : $u^+$(n_q) vs $u^+$ sur le MÊME maillage"
    + ref_str + "\n(élimine l'erreur de maillage, isole l'effet de $n_q$)",
    fontsize=12
)
ax.set_xlim(theta.min(), theta.max())
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9, ncol=2, loc='best', title=r'Ordre de quadrature')

fig.tight_layout()

png_path = os.path.join(save_dir, 'erreur_u_theta_vs_nq_quadrature_pure.png')
fig.savefig(png_path, dpi=200, bbox_inches='tight')
print('saved', png_path)

plt.show()