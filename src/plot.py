import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from abako_utils.read import read_rad


plt.rcParams.update({"font.family": "serif", "font.size": 8})

# Ar main emmisivity lines
HE_ALPHA = (3050, 3160)
LY_ALPHA = (3250, 3350)
HE_BETA = (3590, 3750)
HE_GAMMA = (3800, 4000)
LY_GAMMA = (4100, 4200)

# Default styles
default_style1 = {"linewidth": 1, "color": "k", "label": "full"}
default_style2 = {
    "linewidth": 3,
    "color": mcolors.CSS4_COLORS["dodgerblue"],
    "alpha": 0.7,
    "label": "RDCA",
}

fs_rad_df = read_rad("data/em_op_fs.txt")
rdca_rad_df = read_rad("data/em_op_rdca.txt")


def rad_comparison(
    df1: pd.DataFrame, df2: pd.DataFrame, limits: tuple | None = None, **kwargs
) -> None:
    if not limits:
        limits = (3000, 5000)

    if not "style1" in kwargs:
        style1 = default_style1

    if not "style2" in kwargs:
        style2 = default_style2

    max_fs = df1[df1["E"].between(3590, 3750)]["bb"].max()
    max_rdca = df2[df2["E"].between(3590, 3750)]["bb"].max()

    fs_bb = df1["bb"] #/ max_fs
    rdca_bb = df2["bb"] #/ max_rdca

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6.5), tight_layout=True, dpi=300)

    ax1.plot(df2["E"], rdca_bb, **style2)
    ax1.plot(df1["E"], fs_bb, **style1)
    ax1.legend()

    ax2.plot(df2["E"], rdca_bb, **style2)
    ax2.plot(df1["E"], fs_bb, **style1)
    ax2.set_yscale("log")
    ax2.set_xlabel("Photon energy (eV)")
    ax2.legend()

    ax1.set_xlim(limits)
    ax2.set_xlim(limits)
    fig.supylabel(
        r"Bound-bound emissivity (relative to He-$\beta$)",
        fontsize=plt.rcParams["axes.labelsize"],
    )

    figname = input("Save as: ")
    if not figname:
        figname = "test.png"
    fig.savefig(figname)
