import numpy as np
import pandas as pd

import astropy.units as u
from astropy.constants import R_sun

import paths


def flare_magnetic_field(em, n0, T):
    """Calculte the magnetic field strength of an X-ray flare in Gauss.
    
    Parameters
    ----------
    em : float
        Emission measure in cm^-3.
    n0 : float
        Electron density in cm^-3.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Magnetic field strength in Gauss.
    """

    return 50 * (em / 1e48)**(-0.2) * (n0 / 1e9)**(0.3) * (T / 1e7)**(1.7)

def flare_loop_size(em, n0, T):
    """Calculte the loop size of an X-ray flare in cm.
    
    Parameters
    ----------
    em : float
        Emission measure in cm^-3.
    n0 : float
        Electron density in cm^-3.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Loop size in cm.
    """
    return 1e9 * (em / 1e48)**(0.6) * (n0 / 1e9)**(-0.4) * (T / 1e7)**(-1.6)

def flare_loop_size_from_duration(duration, T, psi=1.2):
    """Calculte the loop size of an X-ray flare in cm.
    
    Parameters
    ----------
    duration : float
        Duration from flare start to peak of emission in ks.
    T : float
        Peak temperature in K.

    Returns
    -------
    float
        Loop size in cm.
    """
    return 0.6 * psi**2 * np.sqrt(T) * duration 

if __name__ == "__main__":

    n0s = np.array([1e11, 1e12, 1e13])

    df = pd.read_csv('../data/mcmc_results.csv', index_col=0)

    # use stellar radius as before
    radius = (0.145 * R_sun).to(u.cm).value

    # calculate B and L
    for n0 in n0s:
        df[f"L_{n0:.0e}"] = flare_loop_size(10**df["EM2_50"], n0, df["T2_50"] * 1e6) / radius
        df[f"B_{n0:.0e}"] = flare_magnetic_field(10**df["EM2_50"], n0, df["T2_50"] *1e6)

    # ------------------------------------------------------------------------------
    # lots of formatting to get the table in the right shape

    # select only the flaring row from df
    flare = df.iloc[2]

    # select only the columns with the loop sizes and B-fields
    flare = flare[[c for c in flare.index if "L_" in c or "B_" in c]]

    flare = flare.reset_index()

    flare["sel"] = flare["index"].apply(lambda x: x.split("_")[0])
    flare["sel2"] = flare["index"].apply(lambda x: x.split("_")[1])

    del flare["index"]

    # add sel as index and sel2 as columns
    flare = flare.pivot(index="sel", columns="sel2")

    # remove the multiindex for both the columns 
    flare.columns = flare.columns.droplevel(0)

    # remove the index name
    flare.index.name = r"$n_0$ [cm$^{-3}$]"

    # prettify
    flare.columns.name = None

    # replace column names with latex
    flare.columns = [r"$10^{" + f"{int(c[-2:])}" + r"}$" for c in flare.columns]


    # rename the index values with latex
    flare.index = [r"$B$ [G]", r"$L$ [R$_\odot$]"]

    # round B row to 0 decimal places
    flare.loc[r"$B$ [G]"] = np.round(flare.loc[r"$B$ [G]"].values.astype(float), 0).astype(int)

    # round L row to 2 decimal places
    flare.loc[r"$L$ [R$_\odot$]"] = np.round(flare.loc[r"$L$ [R$_\odot$]"].values.astype(float), 2).astype(str)

    # print with index name
    string = flare.to_latex(escape=False)

    # layout

    string = string.replace("midrule","hline")
    string = string.replace("toprule","hline")
    string = string.replace("bottomrule","hline")
    string = string.replace("index","")

    print(string)

    # write string to file
    with open(paths.output / "EPIC_flare_loop_table.tex", "w") as f:
        f.write(string)


    # ------------------------------------------------------------------------------
    # Guarcello et al. 2019 method or loop size from rise time
    T = df.loc["flaring", "T2_50"]

    # check out both psi values
    for psi in [1.2, 2.0]:
        print(rf'$\Psi={psi:.1f}$')
        ls = flare_loop_size_from_duration(1e3, 0.13 * (T * 1e6)**1.16, psi=psi) * 100 / radius
        print(fr"${ls:.2f} R_\odot$")