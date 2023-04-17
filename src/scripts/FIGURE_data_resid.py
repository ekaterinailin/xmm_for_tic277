"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in stacked the writefits output from XSPEC, and plots
data, model, and residuals in a two-panel figure.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import paths


def plot_data_resid(file):
    """Plot the data, model and residuals in a two panel figure.
    
    Parameters:
    ----------
    file : str
        filename in results folder
    
    """

    # column names
    names = ["E [keV]", "dE", "flux [counts/s/keV]",
             "e_flux", "model [counts/s/keV]"]

    # read in the file
    _ = pd.read_csv(paths.data / file, delimiter="\s+", skiprows=3, names=names)
    
    # split at the NO NO NO line
    split_at = _[_["E [keV]"]=="NO"].index.values[0]

    # define data and residuals
    data = _.iloc[:split_at].astype(float)
    resid = _.iloc[split_at+1:].drop("model [counts/s/keV]", axis=1).astype(float)

    # make a figure with two subplots, one for data and model, and one for residuals
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7,7),
                           gridspec_kw={'height_ratios': [2.5, 1],
                                        'wspace':0, 'hspace':0})
    # model
    ax[0].stairs(data["model [counts/s/keV]"], 
                 np.append(data["E [keV]"].values-data["dE"].values, 
                           data["E [keV]"].values[-1] + data["dE"].values[-1]), 
                 color="white", edgecolor="k", linewidth=2)

    # flux with errors
    ax[0].errorbar(data["E [keV]"], data["flux [counts/s/keV]"],
                yerr=data["e_flux"], fmt=".", markersize=13, alpha=1, c="olive")

    # residuals
    ax[1].errorbar(resid["E [keV]"], resid["flux [counts/s/keV]"],
                yerr=resid["e_flux"], fmt=".", markersize=13, c="olive")

    # zero line for residuals
    ax[1].axhline(0,c="k")

    # labels
    ax[1].set_xlabel("E [keV]")
    ax[0].set_ylabel("flux [counts/s/keV]")
    ax[1].set_ylabel("residuals")

    # x limits
    for a in ax:
        first = data["E [keV]"].values[0] - data["dE"].values[0]
        last = data["E [keV]"].values[-1] + data["dE"].values[-1]
        a.set_xlim(first, last)

    # layout
    plt.tight_layout()

    # save to file
    filestub = file.split(".")[0]
    plt.savefig(paths.figures / f"{filestub}_data_resid.png", dpi=300)


if __name__ == "__main__":

    # filename
    file = "pn_onlyflare1.txt"

    # plot
    plot_data_resid(file)