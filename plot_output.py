import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.interpolate import griddata

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


def plot_pressure_field(P_values, title_text, suffix):
  # Interpolation des valeurs sur la grille
  Pq = griddata((x, y), P_values, (Xq, Yq), method='linear')

  # Masquer (mettre à NaN) les valeurs à l'intérieur des cercles
  if maillage_data.shape[1] == 3:
    for i in range(maillage_data.shape[0]):
      xc = maillage_data[i, 0]
      yc = maillage_data[i, 1]
      r = maillage_data[i, 2]

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
      xc = maillage_data[i, 0]
      yc = maillage_data[i, 1]
      r = maillage_data[i, 2]

      circle = Circle((xc, yc), r, color='k', fill=True)
      ax.add_patch(circle)
  else:
    print(
        'Le fichier maillage.txt ne contient pas des cercles sous la forme (xc,'
        ' yc, r)'
    )  #

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

