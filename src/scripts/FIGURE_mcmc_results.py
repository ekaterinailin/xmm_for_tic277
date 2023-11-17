"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in the chains from the MCMC analysis, makes corner plots, 
extarct the 16, 50, and 84 percentiles, calculates the weighted mean coronal
T, and plot diagnostic plots, i.e.,

- T1 vs T2
- mean weighted T vs. ratio of emission measures
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import paths
import corner

if __name__ == "__main__":

    # MCMC RESULTS -------------------------------------------------------------

    # read in the MCMC chains and plot them, extract percentiles and save to file

    subsets = [("full data set", "chain_joint_vapec_feo06.fits","x"),
            ("quiescent", "chain_joint_vapec_feo06_noflare.fits","o"),
            ("flaring", "chain_joint_vapec_feo06_flareonly.fits","d"),]

    res = dict()

    for subset, fn, symb in subsets:

        # read in the chain
        df = Table.read(paths.data / fn, format='fits').to_pandas()

        
        # discard the first 5000 steps
        df = df.iloc[5000:]

        norm_to_EM = lambda x: np.log10(1e14 * np.pi * 4. * (13.7 * 3.08567758 * 1e18)**2 * x) #unit cm^-3
        l10 = np.log(10)

        def norm_to_EM_err(norm, norm_err, d, derr):
            # print(norm, norm_err, d, derr, (norm_err**2 / (norm * l10)**2 + 4 * derr**2 / (l10 * d)**2))
            return np.sqrt(norm_err**2 / (norm * l10)**2 + 4 * derr**2 / (l10 * d)**2)

        # convert units
        df["EM1"] = norm_to_EM(df["norm__16"])
        df["EM2"] = norm_to_EM(df["norm__32"])

        df["norm__16"] = df["norm__16"] * 1e6 # now in 10^-6
        df["norm__32"] = df["norm__32"] * 1e6 # now in 10^-6
        
        df["kT__1"] = df["kT__1"] * 11.604525 # now in MK
        df["kT__17"] = df["kT__17"] * 11.604525 # now in MK

        # make a corner plot
        corner.corner(df[["kT__1", "norm__16", "EM1", "kT__17", "norm__32", "EM2"]].values,
                    labels=[r"$T_1$ [MK]", r"$norm_1 \cdot 10^6$", r"$\log_{10} EM_1$ [cm$^{-3}$]",
                            r"$T_2$ [MK]",r"$norm_2  \cdot 10^6$", r"$\log_{10} EM_2$ [cm$^{-3}$]"],
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True,
                    title_kwargs={"fontsize": 12},)
                
        # write to file
        plt.tight_layout()
        figname = "corner" + fn[5:-4] + ".png"
        plt.savefig(paths.figures / figname, dpi=300)

        # caclulate the 0.16, 0.5, 0.84 quantiles and add them to a new dataframe
        q = np.quantile(df[["kT__1","norm__16", "EM1","kT__17","norm__32","EM2"]].values,
                        [0.16, 0.5, 0.84], axis=0).T
        
        # add the quantiles to the results dictionary
        res[subset] = dict(zip(["T1_16","T1_50","T1_84",
                                "norm1_16","norm1_50","norm1_84",
                                "EM1_16","EM1_50","EM1_84",
                                "T2_16","T2_50","T2_84",
                                "norm2_16","norm2_50","norm2_84",
                                "EM2_16","EM2_50","EM2_84",
                                "symb"
                                ],q.flatten()))
        
        # add the symbol to the results dictionary
        res[subset]["symb"] = symb

    # convert to dataframe
    res = pd.DataFrame(res).T

    # --------------------------------------------------------------------------

    # DIAGNOSTICS --------------------------------------------------------------

    # T1 vs T2 -----------------------------------------------------------------
    # plot the two temperatures against each other
    plt.figure(figsize=(6,4.5))
    for label, row in res.iterrows():
        plt.errorbar([row["T1_50"]], [row["T2_50"]],
                    xerr=([row["T1_50"] - row["T1_16"]],
                            [row["T1_84"] - row["T1_50"]]),
                        yerr=([row["T2_50"] - row["T2_16"]],
                            [row["T2_84"] - row["T2_50"]]),
                            markersize=10, fmt=row.symb, label=label, color='grey')
    plt.legend(loc=1, fontsize=12, frameon=False)
    plt.xlabel(r"$T_1$ [MK]", fontsize=12)
    plt.ylabel(r"$T_2$ [MK]",   fontsize=12)

    # set tick label fontsize
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.tight_layout()
    plt.savefig(paths.figures / "T1_vs_T2.png", dpi=300)

    # --------------------------------------------------------------------------


    # ERROR PROPAGATION --------------------------------------------------------

    # calculate the weighted mean of temperatures
    res["weighted_mean_T"] = ((res["T1_50"] * res["norm1_50"] +
                            res["T2_50"] * res["norm2_50"]) / 
                            (res["norm1_50"] + res["norm2_50"]))

    # take the mean of the percentiles as error
    for col in ["T1", "T2", "norm1", "norm2"]:
        # calculate the mean uncertainty from 16 and 84 quantiles
        res[f"{col}_err"] = (res[f"{col}_84"] - res[f"{col}_16"]) / 2

    # error propagate the weighted mean
    et1, et2 = res["T1_err"], res["T2_err"]
    en1, en2 = res["norm1_err"], res["norm2_err"]
    n1, n2 = res["norm1_50"], res["norm2_50"]
    t1, t2 = res["T1_50"], res["T2_50"]


    e2 = ((et1 * n1 / (n1 + n2))**2 +
        (et2 * n2 / (n1 + n2))**2 +
        (en1 * n2 * (t2 - t1) / (n1 + n2)**2)**2 +
        (en2 * n1 * (t1 - t2) / (n1 + n2)**2)**2)

    res["e_weighted_mean_T"] = np.sqrt(e2.values.astype(float))

    # calculate the ratio of EMs
    res["norm_ratio"] = res["norm2_50"] / res["norm1_50"]

    # error propagate the ratio
    e2 = (res["norm_ratio"]**2 * (res["norm1_err"] / res["norm1_50"])**2 +
                                (res["norm2_err"] / res["norm2_50"])**2)

    res["e_norm_ratio"] = np.sqrt(e2.values.astype(float))


    # ERROR propagation from norm to EM with distance uncertanty of 0.11 pc

    # calculate the EM error
    res["e_EM1_50"] = norm_to_EM_err(res["norm1_50"].values.astype(float)*1e-6, res["norm1_err"].values.astype(float)*1e-6, 13.7, 0.11)
    res["e_EM2_50"] = norm_to_EM_err(res["norm2_50"].values.astype(float)*1e-6, res["norm2_err"].values.astype(float)*1e-6, 13.7, 0.11)

    # error propagate the EMs
    

    # --------------------------------------------------------------------------


    # NORM RATIOS VS TEMPERATURE --------------------------------------------------

    # plot the weighted mean temperature vs the emission measure ratio
    plt.figure(figsize=(6.5,4.5))
    for label, row in res.iterrows():
        plt.errorbar([row["weighted_mean_T"]], [row["norm_ratio"]], 
                    xerr=[row["e_weighted_mean_T"]], yerr=[row["e_norm_ratio"]],
                    color='grey',marker=row.symb, label=label, markersize=10)
        
    # labels
    plt.xlabel(r"$EM$-weighted mean coronal temperature $T$ [MK]", fontsize=12)
    plt.ylabel(r"emission measure ratio $EM_{hot}$ / $EM_{cool}$", fontsize=12)

    # layout
    plt.legend(loc=2, frameon=False, fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(0, 5.5)
    plt.tight_layout()

    plt.savefig(paths.figures / "EM_weighted_T_vs_norm_ratio.png", dpi=300)

    # -----------------------------------------------------------------------------

    # save to file
    res.to_csv(paths.data / "mcmc_results.csv", index=True)


