"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script pull the data from Wright et al. 2011 and 2016, and compares with
late M dwarfs in Table 4.
"""


import matplotlib.pyplot as plt
import pandas as pd
import paths


if __name__ == "__main__":

    d_ = pd.read_csv(paths.data / 'wright2016.csv')

    d = d_.sort_values(by='V-K(mag)', ascending=True)
    df = d_[d_["V-K(mag)"] > 5]
    dff = d_[d_["V-K(mag)"] < 5]

    fig, ax = plt.subplots(figsize=(7, 5))

    plt.scatter(dff["Rotation_period(days)"],  dff["log(Lx/Lbol)"], c="silver", 
                alpha=1, s=8,
                label=r"Wright et al. 2011, 2016 ($<$M3.5)")


    plt.scatter(df["Rotation_period(days)"], df["log(Lx/Lbol)"], s=30,
                cmap='viridis', c=df["V-K(mag)"], alpha=0.9, marker='o',
                label=r"Wright et al. 2011, 2016 (M3.5 - M5.5)")

    plt.colorbar(label="V-K [mag]", ax=ax)

    s, ss = 80, 300

    plt.scatter([0.19], [-4.04], s=ss, c='k', marker='*')
    plt.scatter([3.3], [-4], s=ss/1.5, c='k', marker='d') 
    plt.scatter([0.46], [-3.9], s=ss/1.5, c='k', marker='d')
    plt.scatter([0.16], [-3.5], s=ss/1.5, c='k', marker='d')
    plt.scatter([0.61], [-3.1], s=ss/1.5, c='k', marker='d')

    plt.scatter([0.19], [-4.04], s=s, c='r', marker='*', label="TIC 277 (M7)")
    plt.scatter([3.3], [-4], s=s/1.5, c='dimgrey', marker='d', label="TRAPPIST-1 (M8)") 
    plt.scatter([0.46], [-3.9], s=s/1.5, c='cyan', marker='d', label="LHS 248 (M6.5)")
    plt.scatter([0.16], [-3.5], s=s/1.5, c='brown', marker='d', label="NLTT 33370 AB (M7)")
    plt.scatter([0.61], [-3.1], s=s/1.5, c='goldenrod', marker='d', label="LP 412-31 (M8)")

    plt.legend(loc='lower left', fontsize=9, frameon=False)

    plt.xscale('log')
    plt.xlabel("Rotation period [d]", fontsize=14)
    plt.ylabel(r"$\log (L_X\,/\,L_{\rm bol})$", fontsize=14)

    plt.tight_layout()
    plt.savefig(paths.figures / 'lx_lbol.png', dpi=300)