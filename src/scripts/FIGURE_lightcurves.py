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
    xray = pd.read_csv(paths.data / "stacked_xray_lightcurve.csv")
    xray.time = xray.time / 3600. / 24.

    # make the figure

    fix, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
                        gridspec_kw={'height_ratios': [2, 3],
                                            'wspace':0, 'hspace':0})
    ax[0].scatter(om.time, om.rate, label="OM", s=2)
    ax[1].plot(xray.time, xray.normalized_flux + 1., label="PN + MOS1 + MOS2", c="olive")
    
    
    # add inset
    axins = ax[1].inset_axes([0.3, 0.55, 0.4, 0.5])
    l, r = 600, 750
    axins.scatter(om.time[l:r], om.rate[l:r], s=1)
    axins.set_xlim(om.time[l], om.time[r])
    axins.set_ylim(0, 3)
    patch, lines = ax[0].indicate_inset_zoom(axins)
    # [line.set(visible=False) for line in lines]


    # layout
    ax[0].set_ylim(0,3)

    for a in ax:
        a.axhline(1, c="grey")
        a.set_xlim(om.time.min(), xray.time.max())
        a.legend(loc=1, frameon=True)


    ax[0].set_ylabel("normalized flux")
    ax[1].set_ylabel("normalized flux")
    ax[1].set_xlabel("time [days]")

    plt.tight_layout()
    plt.savefig(paths.figures / "lightcurves.png", dpi=300)
