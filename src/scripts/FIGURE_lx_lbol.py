import pandas as pd
import matplotlib.pyplot as plt
import paths

from astropy.table import Table

def tau_wright2018(V, Ks, err=False, eV=None, eKs=None):
    """Convective turnover time from Wright et al. 2018 using
    Eq. 5 from that paper.

    Parameters
    ----------
    V : float
        The V magnitude of the star.
    Ks : float
        The Ks magnitude of the star.
    err : bool, optional
        If True, return the error on the Rossby number.
        The default is False.
    eV : float, optional
        The error on the V magnitude of the star.
        The default is None.
    eKs : float, optional
        The error on the Ks magnitude of the star.
        The default is None.

    Returns 
    -------
    tau : float
        The convective turnover time of the star.
    tau_err_high : float
        The upper error on the convective turnover time.
    tau_err_low : float
        The lower error on the convective turnover time.

    """

    tau = 0.64 + 0.25 * (V - Ks)

    if err:
        tau_err_high = 0.74 + 0.33 * (V + eV - Ks + eKs)
        tau_err_low = 0.54 + 0.17 * (V - eV - Ks - eKs)
        return 10**tau, 10**tau_err_high, 10**tau_err_low
    else:
        return 10**tau



if __name__ == "__main__":

    # read Wright et al. 2011 data
    wright2011 = Table.read(paths.data / "wright2011.fit")
    wright2011["tau"] = tau_wright2018(wright2011["Vmag"], wright2011["Vmag"] - wright2011["V-K"])
    wright2011["Ro"] = wright2011["Prot"] / wright2011["tau"]   

    # read in TIC 277 data
    df = pd.read_csv(paths.data / "apec2_results.csv")

    plt.figure(figsize=(6.5,5))
    
    # TIC 277 data
    df = df[df["label"] == "PN (all)"]
    plt.errorbar(df.Rossby, df.Lx_Lbol, 
                xerr=[df.Rossby-df.Rossby_low,df.Rossby_high-df.Rossby], 
                yerr=df.e_Lx_Lbol, fmt='o', label="TIC 277", alpha=1, c="k")
    
    # Wright et al. 2011 data
    plt.scatter(wright2011["Ro"].value, 10**wright2011["Lx_bol"].value,
                label="Wright et al. 2011", c="grey", zorder=-10, marker="x", alpha=0.5)

    # layout
    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Rossby number", fontsize=14)
    plt.ylabel("L$_x$/L$_{bol}$", fontsize=14)

    plt.legend(fontsize=12, frameon=False, loc=3)
    plt.tight_layout()

    plt.savefig(paths.figures / "lx_lbol.png", dpi=300)
