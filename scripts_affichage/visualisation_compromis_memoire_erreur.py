import numpy as np
import matplotlib.pyplot as plt

# ---------- Mesh h = 0.005 ----------
err_1 = np.array([0.00533617, 7.39839e-05, 9.18779e-06, 1.689e-06,
                  3.15039e-07, 1.57252e-07, 3.28451e-08])
size_H_1 = np.array([9.088242, 9.236115, 9.347534, 9.482101,
                     9.615005, 9.782608, 10.647064])  # MiB
size_dense_1 = 47.211929
tol_1 = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])

# ---------- Mesh h = 0.001 ----------
err_2 = np.array([0.00794499, 0.00224236, 1.07681e-05, 2.0909e-06,
                  2.52421e-07, 3.71211e-08, 8.5019e-09, 2.4422e-09,
                  2.12748e-09, 2.12249e-09, 2.12245e-09, 2.12245e-09])
size_H_2 = np.array([34.372574, 36.375244, 38.123154, 40.071136,
                     41.669159, 43.708466, 47.533951, 57.587646,
                     72.263245, 73.408661, 73.862656, 74.079437])  # MiB
size_dense_2 = 1180.835098
tol_2 = np.array([10.0**(-i) for i in range(2, 14)])  # 12 values

# ---------- Memory saving ----------
saving_1 = 100 * (1 - size_H_1 / size_dense_1)
saving_2 = 100 * (1 - size_H_2 / size_dense_2)

# ---------- Plot ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

datasets = [
    (axes[0], tol_1, saving_1, err_1, "h = 0.005"),
    (axes[1], tol_2, saving_2, err_2, "h = 0.001"),
]

for ax, tol, saving, err, title in datasets:
    # Left axis: memory saving
    l1, = ax.plot(tol, saving, "o-", color="C0", label="taux d'économie mémoire")
    ax.set_xscale("log")
    ax.invert_xaxis()                       # looser tolerance on the left
    ax.set_xlabel("tolérance ACA")
    ax.set_ylabel("Taux d'économie mémoire (%)", color="C0")
    ax.tick_params(axis="y", labelcolor="C0")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.set_title(title)

    # Right axis: relative error
    ax2 = ax.twinx()
    l2, = ax2.plot(tol, err, "s--", color="C3", label="Erreur relative")
    ax2.set_yscale("log")
    ax2.set_ylabel("Erreur relative", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

    ax.legend(handles=[l1, l2], loc="center left")

#fig.suptitle("H-matrix: memory saving and relative error vs. ACA tolerance")
plt.tight_layout()
plt.savefig("memory_saving_error_vs_tolerance.png", dpi=200)
plt.show()

