"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in the flare table and FFD fitting results, and plots the FFD
and power law.
"""


import paths

import numpy as np
import pandas as pd

from altaipony.ffd import FFD

import matplotlib.pyplot as plt

from astropy.modeling import models
from astropy import units as u

if __name__ == "__main__":

    # stellar Teff
    teff = 2680

    # stellar radius
    radius = 0.145

    # TESS band response
    tessresp = pd.read_csv(paths.data / "TESS_response.csv")

    # blackbody
    bb = models.BlackBody(temperature=teff * u.K)

    # blackbody flux in TESS band
    bbtess = bb(tessresp["WAVELENGTH"].values * u.nm)  * tessresp["PASSBAND"].values

    # flux in TESS band
    fstar = np.trapz(bbtess.value, tessresp["WAVELENGTH"].values)

    # flux in all bands
    wav = np.linspace(5,3e4, 20000)
    fall = np.trapz(bb(wav * u.nm).value, wav)

    # luminosity in TESS band
    lumtess = (fstar / fall * bb.bolometric_flux * np.pi * (radius * u.Rsun)**2).to("erg/s")

    # get FFD values
    ffd_vals = pd.read_csv(paths.data / "tess_ffd.csv")

    # get flares
    df = pd.read_csv(paths.data / "tess_flares.csv")

    # convert ED to E
    df["ed_rec"] = df["ed_rec"] * lumtess
    df["ed_rec_err"] = df["ed_rec_err"] * lumtess

    # make FFD
    ffd = FFD(df, tot_obs_time=df["tot_obs_time"].values[0])
    ffd.alpha = ffd_vals["alpha"].values[0]
    ffd.beta = ffd_vals["beta"].values[0] * lumtess**(ffd.alpha -1)
    ffd.alpha_low_err = ffd_vals["alpha_low_err"].values[0]
    ffd.alpha_up_err = ffd_vals["alpha_up_err"].values[0]
    ed, freq, counts = ffd.ed_and_freq()

    # plot
    fig, ax = plt.subplots(1,1,figsize=(5,4))

    # make label
    label = (fr"$\alpha =$ {ffd.alpha:.1f} (+{ffd.alpha_up_err:.1f} "
            fr"/ -{ffd.alpha_low_err:.1f})")

    # FFD
    plt.scatter(ed, freq, c="k")

    # power law
    ffd.plot_powerlaw(ax, c="olive", label=label)

    # layout
    plt.legend(loc=1, frameon=False)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$E_{flare}$ [erg]")
    plt.ylabel(r"flares per day above $E_{flare}$")
    plt.tight_layout()

    # save to file
    plt.savefig(paths.figures / "tess_ffd.png", dpi=300)