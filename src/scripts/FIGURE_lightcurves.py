"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in stacked PN and MOS, and optical monitoring light curves 
from XMM-Newton, and plots them together.
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import paths


if __name__ == "__main__":

    # read in optical data
    om = pd.read_csv(paths.data / "timeseries.csv")
    om.time = om.time / 3600. / 24.
    om.rate = om.rate / np.nanmedian(om.rate)

    # read in X-ray data
    xray = pd.read_csv(paths.data / "corrected_merged_epic_lc.csv")
    xray["TIME"] = xray["TIME"] / 3600. / 24.

    # make the figure

    fix, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
                        gridspec_kw={'height_ratios': [2, 3],
                                            'wspace':0, 'hspace':0})
    ax[0].scatter(om.time, om.rate, label="OM", s=2, alpha=0.8)
    ax[1].errorbar(xray["TIME"], xray["RATE"], yerr=xray["ERROR"],
                   label="PN + MOS1 + MOS2", c="olive", alpha=0.8, lw=1.5)
    
    # add vertical filled area
    ax[1].axvspan(776075500 / 3600. / 24., 776080000 / 3600. / 24.,
                  color="grey", alpha=0.2)
    
    # add inset
    axins = ax[1].inset_axes([0.3, 0.55, 0.4, 0.5])
    l, r = 600, 750
    axins.scatter(om.time[l:r], om.rate[l:r], s=1)
    axins.set_xlim(om.time[l], om.time[r])
    axins.set_ylim(0, 3)
    patch, lines = ax[0].indicate_inset_zoom(axins)
   

    # layout
    ax[0].set_ylim(0,3)
    ax[0].axhline(1, c="grey")

    for a in ax:
        
        a.set_xlim(om.time.min(), xray["TIME"].max())
        a.legend(loc=1, frameon=True)


    ax[0].set_ylabel("normalized flux")
    ax[1].set_ylabel("background subtracted flux <10 keV [cts/s]")
    ax[1].set_xlabel("time [days]")

    plt.tight_layout()
    plt.savefig(paths.figures / "lightcurves.png", dpi=300)
