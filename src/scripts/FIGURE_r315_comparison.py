"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script plots the FFD results from Medina et al. 2020, Murray et al. 2022,
Trappist-1, and the four stars with high latitude flares from Ilin et al. 2021.
"""

import paths

import pandas as pd
import numpy as np

from astropy.table import Table

import matplotlib.pyplot as plt

import adjustText as aT

if __name__ == "__main__":

    # get Medina et al table
    df = Table.read(paths.data / "medina2020.fit")
    df = df.to_pandas()

    # get FFD values
    ffd_vals = pd.read_csv(paths.data / "tess_ffd.csv")

    # get rotation periods
    ilin2021 = pd.read_csv(paths.data / "ilin2021updated_w_Rossby_Lbol.csv")

    ffd_vals = ffd_vals.merge(ilin2021, on="TIC")

    # calculate r315
    r315s = []
    for tic, ff in ffd_vals.groupby("TIC"):

        # energy beta
        with open(paths.data / f"{tic}_energy_beta.txt", "r") as f:
            beta = float(f.read())

        f315 = beta / (ff.alpha - 1) * (10**31.5)**(-ff.alpha +1)
        r315s.append(np.log10(f315).values[0])

    ffd_vals["r315"] = r315s


    # TRAPPIST-1
    trapp_r315 = np.log10(10**(-2.5) * 24) # paudel 2018
    trapp_ro = 0.047 # roettenbacher 2017
    trapp_rot = 3.3 

    # data from Murray 2022 Table 1
    murray = [("M4-M5",np.log10(.055),"dashed"),
            ("M6", np.log10(0.112),"dotted"),
            ("M7", np.log10(0.06),"solid")]


    # make plot
    plt.figure(figsize=(8,4))

    # add Murray+2022
    for label, r, sty in murray:
        plt.axhline(r, label=label + " (Murray+2022)", linestyle=sty, c="darkgrey",zorder=-10)


    # color code Medina+2020 by mass
    color_marker_label = [("palegreen","o",0.1,0.12),
                        ("yellow","d",0.12,0.15),
                        ("green","*",0.15,0.2),
                        ("cyan","x",0.2,0.25)]

    # plot Medina+2020
    for c, m, low, high in color_marker_label:
        g = df.loc[(df.Mstar>low) & (df.Mstar<high)]
        plt.scatter(g.Prot, g.Rate, c="k", marker=m, alpha=1, s=45)
        plt.scatter(g.Prot, g.Rate, c=c, marker=m, alpha=1,s=30,
                    label=rf"${low:.2f}<M_*<{high:.2f}M_\odot$ (Medina+2020)")

    # plot TRAPPIST-1
    plt.scatter([trapp_rot],[trapp_r315],marker=rf"$T$", s=60, 
                label="TRAPPIST-1 (Paudel+2018)", c="k")

    # plot our values and annotate with SpT
    txts = []
    for i, d in ffd_vals.iterrows():

        plt.scatter([d.Prot_days],[d.r315],c="k",s=35,marker="s")
        plt.scatter([d.Prot_days],[d.r315],c="magenta",s=20,marker="s")
        
        txts.append(plt.text( d.Prot_days, d.r315, d.SpT,
                        fontsize=11,  color="k"))
        
    # legend and layout    
    plt.legend(loc=(1.05,0.3), frameon=False, )

    plt.xscale("log")

    aT.adjust_text(txts, )
    plt.xlabel("rotation period [d]")
    plt.ylabel(r"log$_{10}$ flares per day above log$_{10}\,E = 31.5$ erg")
    plt.ylim(-14,1)

    plt.tight_layout()

    # save
    plt.savefig(paths.figures / "r315_prot.png", dpi=300)