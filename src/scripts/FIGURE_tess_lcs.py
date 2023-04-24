"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in the de-trended TESS light curves and plots them in a
three-panel figure.
"""


from astropy.table import Table
import matplotlib.pyplot as plt
import paths

import numpy as np


if __name__ == "__main__":

    # define sectors
    sectors = 12, 37, 39

    # read in light curves
    lcrs = [Table.read(paths.data / f"tic277_tess_detrended_{s}.fits") for s in sectors]

    # make three subplots with two light curves each, one showing the flux the other the
    # detrended flux
    fig, axes = plt.subplots(3, 1, figsize=(10, 9))

    # loop over axes, light curves and sectors
    for ax, lcr, sector in zip(axes, lcrs, sectors):

        # plot the un-detrended and detrended flux
        ax.plot(lcr['TIME'], lcr['FLUX']/np.nanmedian(lcr["FLUX"]), 'k.',
                label=f"sector {sector}", ms=1.5)
        ax.plot(lcr['TIME'], lcr['DETRENDED_FLUX']/np.nanmedian(lcr["DETRENDED_FLUX"]) +0.2, 
                'r.', ms=1.5)
        
        # layout
        ax.set_ylabel("flux")
        ax.set_ylabel("normalized flux")
        ax.set_ylim(0.9, 1.5)
        ax.set_xlim(lcr["TIME"][0], lcr["TIME"][-1])

        # add sector number
        ax.text(0.02, 0.9, f"Sector {sector}", transform=ax.transAxes)

    # shared x-axis
    axes[-1].set_xlabel("time [BJD - 2457000]")

    plt.tight_layout()
    plt.savefig(paths.figures / "tess_lcs.png", dpi=250)