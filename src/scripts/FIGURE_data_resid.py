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


def plot_data_resid_3(file):


     # column names
    names = ["E [keV]", "dE", "flux [counts/s/keV]",
             "e_flux", "model [counts/s/keV]"]

    # read in the file
    _ = pd.read_csv(paths.data / file, delimiter="\s+", skiprows=3, names=names)
    
    # split at the NO NO NO line
    split_at = _[_["E [keV]"]=="NO"].index.values

    d = _[:split_at[0]], _[split_at[0]+1:split_at[1]], _[split_at[1]+1:split_at[2]]
    r = _[split_at[2]+1:split_at[3]], _[split_at[3]+1:split_at[4]], _[split_at[4]+1:]

    # make a figure with two subplots, one for data and model, and one for residuals
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7,7),
                           gridspec_kw={'height_ratios': [2.5, 1],
                                        'wspace':0, 'hspace':0})
    
    for din, rin, c, label, m in zip(d, r, ["olive", "orange", "blue"],
                                           ["PN", "MOS1", "MOS2"],
                                           ["o", "s", "x"]):

        data = din.astype(float)
        resid = rin.drop("model [counts/s/keV]", axis=1).astype(float)

         # model
        ax[0].stairs(data["model [counts/s/keV]"], 
                    np.append(data["E [keV]"].values-data["dE"].values, 
                            data["E [keV]"].values[-1] + data["dE"].values[-1]), 
                    color="white", edgecolor=c, linewidth=2, label=label)

        # flux with errors
        ax[0].errorbar(data["E [keV]"], data["flux [counts/s/keV]"],
                    yerr=data["e_flux"]/1.644854, fmt=m, markersize=6, alpha=1, c=c)

        # residuals
        ax[1].errorbar(resid["E [keV]"], resid["flux [counts/s/keV]"],
                    yerr=resid["e_flux"]/1.644854, fmt=m, markersize=6, c=c)
        
        
      # zero line for residuals
    ax[1].axhline(0,c="k")

    # labels
        # labels
    ax[1].set_xlabel("E [keV]", fontsize=14)
    ax[0].set_ylabel(r"flux [counts s$^{-1}\,$keV$^{-1}$]", fontsize=13)
    ax[1].set_ylabel(r"residuals [counts s$^{-1}\,$keV$^{-1}$]", fontsize=13)

    # increase tick size labels
    for a in ax:
        a.tick_params(axis='both', which='major', labelsize=13)
    

    # y limits
    ax[0].set_ylim(0,)

    # legend
    ax[0].legend(loc=1, frameon=False, fontsize=12)

    # x limits
    for a in ax:
        # first = data["E [keV]"].values[0] - data["dE"].values[0]
        # last = data["E [keV]"].values[-1] + data["dE"].values[-1]
        a.set_xlim(0.25, 5)

        a.set_xscale("log")

    ax[1].set_xticks(ticks=[0.25,0.4,0.6,0.8,1,2,3,4,5],
                     labels=[0.25,0.4,0.6,0.8,1,2,3,4,5])

    # layout
    plt.tight_layout()

    # save to file
    filestub = file.split(".")[0]
    plt.savefig(paths.figures / f"{filestub}_data_resid.png", dpi=300)



def plot_data_resid_1(file):
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
                yerr=data["e_flux"] / 1.644854, fmt=".", markersize=13, alpha=1, c="olive")

    # residuals
    ax[1].errorbar(resid["E [keV]"], resid["flux [counts/s/keV]"],
                yerr=resid["e_flux"]/1.644854, fmt=".", markersize=13, c="olive")

    # zero line for residuals
    ax[1].axhline(0,c="k")

    # labels
    ax[1].set_xlabel("E [keV]", fontsize=15)
    ax[0].set_ylabel(r"flux [counts s$^{-1}$ keV^$^{-1}$]", fontsize=15)
    ax[1].set_ylabel(r"residuals [counts s$^{-1}$ keV^$^{-1}$]", fontsize=15)

    # increase tick size labels
    ax[0].tick_params(axis='both', which='major', labelsize=15)

    # x limits
    for a in ax:
        first = data["E [keV]"].values[0] - data["dE"].values[0]
        last = data["E [keV]"].values[-1] + data["dE"].values[-1]
        a.set_xlim(first, last)
        a.set_xscale("log")

    # layout
    plt.tight_layout()

    # save to file
    filestub = file.split(".")[0]
    plt.savefig(paths.figures / f"{filestub}_data_resid.png", dpi=300)


if __name__ == "__main__":

    # filename
    # ["pn_onlyflare.txt", "joint_onlyflare.txt", "joint_noflare.txt",
    #          "joint_all.txt", 'mos1.txt', 'mos2.txt', 'pn.txt', "pn_noflare.txt",
    #          "joint_3vapec_fe06.txt", "pn_all_vapec_fe06.txt", 
    #          "pn_onlyflare_3vapec_fe06.txt", "joint_vapec_fe06.txt", 
    #          "pn_noflare_vapec_fe06.txt", "pn_onlyflare_vapec_fe06.txt",
    files =  ["joint_chain_fit.txt"]

    # plot
    for file in files:
        if file.split("_")[0] == "joint":
            plot_data_resid_3(file)
        else:
            plot_data_resid_1(file)