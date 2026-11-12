import numpy as np
import matplotlib.pyplot as plt

def plot_gain_error(k, gain, err, title, fname, ref=None):
    fig, ax = plt.subplots(figsize=(7, 5))

    l1, = ax.plot(k, gain, "o-", color="C0", label="Memory gain")
    handles = [l1]

    # Optional reference curve (dotted grey)
    if ref is not None:
        lref, = ax.plot(k, ref, ":", color="grey", linewidth=2,
                        label="$1 - h\\,|\\log h|$")
        handles.append(lref)

    ax.set_xscale("log")
    ax.set_xlabel("Wavenumber $k$")
    ax.set_ylabel("Memory gain  $1 - \\mathrm{size}_{\\mathcal{H}}/\\mathrm{size}_{dense}$  (%)", color="C0")
    ax.tick_params(axis="y", labelcolor="C0")
    ax.grid(True, which="both", ls=":", alpha=0.6)

    ax2 = ax.twinx()
    l2, = ax2.plot(k, err, "s--", color="C3", label="Relative error")
    ax2.set_yscale("log")
    ax2.set_ylabel("Relative error", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")
    handles.append(l2)

    ax.legend(handles=handles, loc="center left")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    plt.show()

# ---------- Dataset 1 (unchanged) ----------
mem_H_1 = np.array([53.344803, 72.263245, 73.643265, 73.744568, 73.996521,
                    74.238586, 74.231232, 74.342804, 75.647629, 75.647629])
mem_dense_1 = 1180.835098
err_1 = np.array([2.43785e-10, 2.12748e-09, 7.60085e-06, 0.371511, 0.559874,
                  0.204963, 0.0297106, 0.00872098, 0.00120242, 0.000419408])
k_1 = np.array([10.0**i for i in range(0, 10)])
gain_1 = 100 * (1 - mem_H_1 / mem_dense_1)
plot_gain_error(k_1, gain_1, err_1,
                "h = 0.001, ACA tolerance = 1e-10",
                "gain_error_vs_k_fixed_h.png")

# ---------- Dataset 2: kh = 10 constant, with reference curve ----------
k_2 = np.array([10, 100, 200, 500, 1000, 2000])
h_2 = np.array([0.1, 0.01, 0.005, 0.002, 0.001, 0.0005])
mem_H_2 = np.array([0.118164, 3.602966, 10.778397, 31.814667, 73.744568, 167.188370])
mem_dense_2 = np.array([0.118164, 11.816406, 47.211929, 295.275894, 1180.835098, 4723.877335])
err_2 = np.array([8.83746e-17, 1.74314e-05, 0.0434265, 0.264187, 0.371511, 0.471906])

gain_2 = 100 * (1 - mem_H_2 / mem_dense_2)
ref_2 = 100 * (1 - h_2 * np.abs(np.log(h_2)))   # = 100 * (1 + h ln h)

plot_gain_error(k_2, gain_2, err_2,
                "$kh = 10$ constant",
                "gain_error_vs_k_fixed_kh.png",
                ref=ref_2)