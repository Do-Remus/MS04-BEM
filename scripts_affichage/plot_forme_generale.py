import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.path import Path
from scipy.interpolate import griddata


# --- Fichiers (produits par le bloc TP_balle_golf de main.cpp) ---
data_file = '../outputs/narval_u_k10.000000_R1.000000_L10.000000_hM0.000500_hS0.050000_d0.000005_N25_q4.txt'
mesh_file = '../outputs/narval_maillage.txt'

# --- Champ à tracer : 'diffracte' (u+) ou 'total' (u+ + u_inc) ---
FIELD = 'diffracte'

# --- Paramètres d'affichage du maillage ---
MESH_COLOR = 'blue'
MESH_LINEWIDTH = 1.0
NORMAL_COLOR = 'red'
NORMAL_STEP = 10        # une normale toutes les NORMAL_STEP segments
NORMAL_LENGTH = 0.15    # longueur des flèches (unités du graphe)
GRID_N = 300            # résolution de la grille d'interpolation

# 1. Données : x y Re(u) Im(u) Re(u+uinc) Im(u+uinc)
data = np.loadtxt(data_file)
x, y = data[:, 0], data[:, 1]
if FIELD == 'total':
    P_real, P_imag = data[:, 4], data[:, 5]
    field_name = 'u+ + u_inc'
else:
    P_real, P_imag = data[:, 2], data[:, 3]
    field_name = 'u+'

base_name = os.path.splitext(os.path.basename(data_file))[0] + '_' + FIELD
save_dir = '../outputs'
os.makedirs(save_dir, exist_ok=True)

# Grille régulière pour l'interpolation
xq = np.linspace(x.min(), x.max(), GRID_N)
yq = np.linspace(y.min(), y.max(), GRID_N)
Xq, Yq = np.meshgrid(xq, yq)

# 2. Maillage (x1 y1 x2 y2 par ligne)
mesh = np.atleast_2d(np.loadtxt(mesh_file))
x1, y1, x2, y2 = mesh[:, 0], mesh[:, 1], mesh[:, 2], mesh[:, 3]

# Polygone de la balle (sommets P1 des segments, dans l'ordre)
ball_poly = np.column_stack([x1, y1])
ball_path = Path(ball_poly)
inside_ball = ball_path.contains_points(
    np.column_stack([Xq.ravel(), Yq.ravel()])).reshape(Xq.shape)


def plot_mesh(ax):
    """Superpose les segments du maillage et leurs normales sur `ax`."""
    segments = np.stack([np.column_stack([x1, y1]),
                         np.column_stack([x2, y2])], axis=1)
    ax.add_collection(LineCollection(segments, colors=MESH_COLOR,
                                     linewidths=MESH_LINEWIDTH, zorder=3))

    idx = np.arange(0, len(x1), NORMAL_STEP)
    dx = x2[idx] - x1[idx]
    dy = y2[idx] - y1[idx]
    norm = np.hypot(dx, dy)
    ok = norm > 0
    idx, dx, dy, norm = idx[ok], dx[ok], dy[ok], norm[ok]

    mid_x = x1[idx] + dx / 2
    mid_y = y1[idx] + dy / 2
    u = NORMAL_LENGTH * dy / norm     # normale sortante (balle orientée anti-horaire)
    v = -NORMAL_LENGTH * dx / norm
    #ax.quiver(mid_x, mid_y, u, v, angles='xy', scale_units='xy', scale=1,
     #         color=NORMAL_COLOR, width=0.004, zorder=5)


def plot_pressure_field(P_values, title_text, suffix):
    Pq = griddata((x, y), P_values, (Xq, Yq), method='linear')
    Pq[inside_ball] = np.nan   # masque l'intérieur de la balle

    fig, ax = plt.subplots(figsize=(8, 8))
    extent = [xq.min(), xq.max(), yq.min(), yq.max()]
    im = ax.imshow(Pq, extent=extent, origin='lower', aspect='equal',
                   cmap='viridis')

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_title(title_text, fontsize=20)
    ax.tick_params(labelsize=20)

    # Balle pleine (noire) à la place des cercles
    ax.add_patch(Polygon(ball_poly, closed=True, facecolor='k',
                         edgecolor='none', zorder=2))

    plot_mesh(ax)

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    png_path = os.path.join(save_dir, f'{base_name}_{suffix}.png')
    fig.savefig(png_path, dpi=200, bbox_inches='tight')
    print('saved', png_path)
    plt.show()


plot_pressure_field(P_real, f'Partie réelle de {field_name}', 'real')
plot_pressure_field(P_imag, f'Partie imaginaire de {field_name}', 'imag')
plot_pressure_field(np.hypot(P_real, P_imag), f'Norme de {field_name}', 'norm')