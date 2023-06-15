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
from astropy.constants import sigma_sb


def flare_factor(teff, radius, wav, resp,  tflare=10000):
    """Calculate the flare energy factor in ergs.

    Parameters
    ----------
    teff : float
        Stellar effective temperature in Kelvin.
    radius : float
        Stellar radius in solar radii.
    wav : array
        Array of wavelengths in nanometers.
    resp : array
        Array of bandpass responses.
     tflare : float
        Flare temperature in Kelvin.
    
    Returns
    -------
    factor : float
        Flare energy factor in ergs/s.
    """

    # blackbody
    bb = models.BlackBody(temperature=teff * u.K)

    # blackbody flux in TESS band
    bbwavs = bb(wav * u.nm)  * resp

    fluxs = np.trapz(bbwavs.value, wav)

    # blackbody
    bb = models.BlackBody(temperature=tflare * u.K)

    # blackbody flux in TESS band
    bbwavf = bb(wav * u.nm)  * resp

    fluxf = np.trapz(bbwavf.value, wav)

    ratio = fluxs / fluxf

    factor = ratio * np.pi * (radius * u.R_sun) ** 2 * sigma_sb * (tflare * u.K)**4

    return factor.to("erg/s")


if __name__ == "__main__":

    tic_teff_rad = [(277539431, 2680, 0.145),
                    (237880881, 3060, 0.275),
                    (452922110, 2680, 0.137),
                    (44984200, 2810, 0.145)]
    
    

    for tic, teff, radius in tic_teff_rad:

        # TESS band response
        tessresp = pd.read_csv(paths.data / "TESS_response.csv")

        factor = flare_factor(teff, radius, tessresp["WAVELENGTH"].values, tessresp["PASSBAND"].values)

        # get FFD values
        ffd_vals = pd.read_csv(paths.data / "tess_ffd.csv")
        ffd_vals = ffd_vals[ffd_vals["TIC"] == tic]

        # get flares
        df = pd.read_csv(paths.data / "tess_flares.csv")
        df = df[df["TIC"] == tic]

        # convert ED to E
        df["ed_rec"] = df["ed_rec"] * factor
        df["ed_rec_err"] = df["ed_rec_err"] * factor

        # make FFD
        ffd = FFD(df, tot_obs_time=df["tot_obs_time"].values[0])
        ffd.alpha = ffd_vals["alpha"].values[0]
        ffd.beta = ffd_vals["beta"].values[0] * factor**(ffd.alpha -1)
        ffd.alpha_low_err = ffd_vals["alpha_low_err"].values[0]
        ffd.alpha_up_err = ffd_vals["alpha_up_err"].values[0]
        ed, freq, counts = ffd.ed_and_freq()

        # plot
        fig, ax = plt.subplots(1,1,figsize=(5,4))

        # make label
        label = (fr"$\alpha =$ {ffd.alpha:.1f} (+{ffd.alpha_up_err:.1f} "
                fr"/ -{ffd.alpha_low_err:.1f})")
        
        # if TIC 277, add OM flare 
        if tic == 277539431:
            om = pd.read_csv(paths.data / "flare_energies.csv")
            om = om[om.instrument == "OM"].iloc[0]
            plt.errorbar([om.E_erg], [om.rate_per_day], 
                         yerr = [[0.5*om.rate_per_day], [om.rate_per_day]],
                          xerr=om.eE_erg, c="grey", marker="s", label="OM") 

            x = np.linspace(om.E_erg/10,np.min(ed),10)
            pl = ffd.beta / (ffd.alpha - 1) * x**(- ffd.alpha + 1)

            plt.plot(x, pl, linestyle="dashed", c= "grey")

        # FFD
        plt.scatter(ed, freq, c="k", label="TESS")

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
        plt.savefig(paths.figures / f"{tic}_tess_ffd.png", dpi=300)

        # save new_beta
        with open(paths.data / f"{tic}_energy_beta.txt", "w") as f:
            f.write(str(ffd.beta.value))