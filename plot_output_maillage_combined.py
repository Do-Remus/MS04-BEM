import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from scipy.interpolate import griddata

# --- Paramètres d'affichage du maillage ---
MESH_COLOR = 'blue'
MESH_LINEWIDTH = 1.5
NORMAL_COLOR = 'red'
NORMAL_STEP = 100       # une normale toutes les NORMAL_STEP segments
NORMAL_LENGTH = 1.0    # longueur des flèches (unités du graphe)

# 1. Charger les données depuis le fichier
data_file = 'outputs/u_N_plusResu_k=10.000000_a=1.000000_N=15.txt'
data = np.loadtxt(data_file)

# Nom de base des PNG : même nom que le fichier de données (contient k, a et N)
base_name = os.path.splitext(os.path.basename(data_file))[0]
save_dir = 'outputs'
os.makedirs(save_dir, exist_ok=True)

# Extraire les colonnes
x = data[:, 0]
y = data[:, 1]
P_real = data[:, 2]
P_imag = data[:, 3]

# Créer une grille régulière pour interpolation
xq = np.linspace(x.min(), x.max(), 100)
yq = np.linspace(y.min(), y.max(), 100)
Xq, Yq = np.meshgrid(xq, yq)

# Charger les obstacles et s'assurer qu'ils sont en 2 dimensions
maillage_data = np.loadtxt('outputs/obstacles.txt')
maillage_data = np.atleast_2d(maillage_data)

# Charger le maillage (segments) une seule fois
mesh = np.atleast_2d(np.loadtxt('outputs/test.txt'))
x1, y1, x2, y2 = mesh[:, 0], mesh[:, 1], mesh[:, 2], mesh[:, 3]


def plot_mesh(ax):
    """Superpose les segments du maillage et leurs normales sur `ax`."""
    # Segments (tracés d'un coup avec un LineCollection, beaucoup plus rapide)
    segments = np.stack([np.column_stack([x1, y1]),
                         np.column_stack([x2, y2])], axis=1)
    ax.add_collection(LineCollection(segments, colors=MESH_COLOR,
                                     linewidths=MESH_LINEWIDTH, zorder=3))
    # Points aux extrémités
    ax.plot(np.concatenate([x1, x2]), np.concatenate([y1, y2]), linestyle='none',
            marker='o', markersize=2, markerfacecolor='black',
            markeredgecolor=MESH_COLOR, zorder=4)

    # Normales (sous-ensemble des segments)
    idx = np.arange(0, len(x1), NORMAL_STEP)
    dx = x2[idx] - x1[idx]
    dy = y2[idx] - y1[idx]
    norm = np.hypot(dx, dy)
    ok = norm > 0  # évite la division par zéro
    idx, dx, dy, norm = idx[ok], dx[ok], dy[ok], norm[ok]

    mid_x = x1[idx] + dx / 2
    mid_y = y1[idx] + dy / 2
    u = NORMAL_LENGTH * dy / norm
    v = -NORMAL_LENGTH * dx / norm
    ax.quiver(mid_x, mid_y, u, v, angles='xy', scale_units='xy', scale=1,
              color=NORMAL_COLOR, width=0.005, zorder=5)


def plot_pressure_field(P_values, title_text, suffix):
    # Interpolation des valeurs sur la grille
    Pq = griddata((x, y), P_values, (Xq, Yq), method='linear')

    # Masquer (mettre à NaN) les valeurs à l'intérieur des cercles
    if maillage_data.shape[1] == 3:
        for i in range(maillage_data.shape[0]):
            xc, yc, r = maillage_data[i, :3]
            inside_circle = (Xq - xc) ** 2 + (Yq - yc) ** 2 <= r**2
            Pq[inside_circle] = np.nan

    fig, ax = plt.subplots(figsize=(8, 8))

    # Affichage de la carte de chaleur
    extent = [xq.min(), xq.max(), yq.min(), yq.max()]
    im = ax.imshow(
        Pq, extent=extent, origin='lower', aspect='equal', cmap='viridis'
    )

    # Colorbar et axes
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_title(title_text, fontsize=20)
    ax.tick_params(labelsize=20)

    # Superposition des cercles noirs pleins
    if maillage_data.shape[1] == 3:
        for i in range(maillage_data.shape[0]):
            xc, yc, r = maillage_data[i, :3]
            ax.add_patch(Circle((xc, yc), r, color='k', fill=True))
    else:
        print(
            'Le fichier obstacles.txt ne contient pas des cercles sous la forme'
            ' (xc, yc, r)'
        )

    # Superposition du maillage et des normales
    plot_mesh(ax)

    # Garder les limites de la carte de chaleur (le maillage ne les modifie pas)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    # Sauvegarde du PNG (avant plt.show(), sinon la figure est vidée à la fermeture)
    png_path = os.path.join(save_dir, f'{base_name}_{suffix}.png')
    fig.savefig(png_path, dpi=200, bbox_inches='tight')
    print('saved', png_path)

    plt.show()


# --- 1. Partie réelle de la pression P ---
plot_pressure_field(P_real, 'Partie réelle de la solution exterieure u+', 'real')

# --- 2. Partie imaginaire de la pression P ---
plot_pressure_field(P_imag, 'Partie imaginaire de la solution exterieure u+', 'imag')

# --- 3. Norme de la pression complexe P ---
P_norm = np.sqrt(P_real**2 + P_imag**2)
plot_pressure_field(P_norm, 'Norme de la solution exterieure complexe u+', 'norm')