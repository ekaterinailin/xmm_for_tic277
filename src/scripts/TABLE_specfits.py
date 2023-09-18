"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in the XPEC fit results and derived parameters and produces
a LaTeX table.
"""

import pandas as pd
import numpy as np
import paths

def convert_to_scinote(series, rel_err=1e-2):
    """Convert a series to scientific notation.

    Parameters
    ----------
    series : pd.Series
        The series to convert.

    Returns
    -------
    pd.Series
        The converted series of strings.
    """
    return series.apply(lambda x: f"{x:.1e}" if np.abs(x)<rel_err else f"{x:.2f}").astype(str)



def tex_up_low(val, err, err2):
    """Convert a value and two errors into a LaTeX string.

    Parameters
    ----------
    val : float
        The value.
    err : float
        The first error.
    err2 : float
        The second error.
    
    Returns
    -------
    str
        The LaTeX string.
    """

    return ("$" +
            convert_to_scinote(val) +
            "^{" + 
            convert_to_scinote(err) + 
            "}_{" +
            convert_to_scinote(err2) + "}$")

def tex_one_err(val, err):
    """Convert a value and one error into a LaTeX string.

    Parameters
    ----------
    val : float
        The value.
    err : float
        The error.

    Returns
    -------
    str
        The LaTeX string.
    """

    return ("$" +
            convert_to_scinote(val) +
            "[" +
            convert_to_scinote(err) +
            "]$")

if __name__ == "__main__":

    df = pd.read_csv(paths.data / "mcmc_results.csv", index_col=0)

    # load luminosity from PN data
    lxs = pd.read_csv(paths.data / "joint_vapec_chain_fits.csv")

    # get flux and luminosity
    g = lambda x: fr"${x.flux_erg_s_cm2/1e-14:.2e} [{x.flux_erg_s_cm2_err/1e-14:.2e}]$"
    lxs[r"$F_X$ [$10^{-14}$ erg s$^{-1}$ cm$^{-2}$]"] = lxs.apply(g, axis=1)

    g = lambda x: fr"${x.Lx_erg_s/1e26:.2f} [{x.Lx_erg_s_err/1e26:.2f}]$"
    lxs[r"$L_X$ [$10^{26}$ erg s$^{-1}$]"] = lxs.apply(g, axis=1) 

    # rename subset in lxs as follows: spec to full data set, 
    # no flare to quiescent, onlyflare to flaring
    lxs["subset"] = lxs["subset"].replace({"spec": "full data set", 
                                           "noflare": "quiescent", 
                                           "onlyflare": "flaring"})

    # merge the subset, and LX column in lxs with df
    df = df.reset_index()
    df = df.merge(lxs[["subset",
                       r"$L_X$ [$10^{26}$ erg s$^{-1}$]",
                       r"$F_X$ [$10^{-14}$ erg s$^{-1}$ cm$^{-2}$]"]], 
                       right_on="subset", left_on="index")
    
    df = df.set_index("index")

    # remove extra columns
    del df["symb"]
    del df["norm_ratio"]
    del df["e_norm_ratio"]
    del df["subset"]

    for c in df.columns:
        if "_err" in c:
            del df[c]

    # rename columns
    mapcolnames = {'T1': r'$T_{\rm cool}$ [MK]',
                   'T2': r'$T_{\rm hot}$ [MK]',
                   'norm1': r'$10^6$ norm$_{\rm cool}$',
                   'norm2': r'$10^6$ norm$_{\rm hot}$',
                   'weighted_mean_T': r'$T_{\rm mean}$ [MK]'}
    
    # convert to LaTeX
    for col in ["T1", "T2", "norm1", "norm2"]:

        df[f"{col}_uperr"] = df[f"{col}_84"] - df[f"{col}_50"]
        df[f"{col}_lowerr"] = df[f"{col}_16"] - df[f"{col}_50"]

        newname = mapcolnames[col]
        df[newname] = tex_up_low(df[f"{col}_50"], 
                                 df[f"{col}_uperr"], 
                                 df[f"{col}_lowerr"])

        del df[f"{col}_84"]
        del df[f"{col}_50"]
        del df[f"{col}_16"]
        del df[f"{col}_uperr"]
        del df[f"{col}_lowerr"]


    col = "weighted_mean_T"
    df[mapcolnames[col]] = tex_one_err(df[col], df[f"e_{col}"])
    del df[col]
    del df[f"e_{col}"]

    string = df.T.to_latex(escape=False,index=True)

    # layout
    string = string.replace("midrule","hline")
    string = string.replace("toprule","hline")
    string = string.replace("bottomrule","hline")
    string = string.replace("index","")

    print(string)
   

    # write to file
    path = paths.output / f"mcmc_specfit.tex"
    print("Write to file: ", path)

    with open(path, "w") as f:
        f.write(string)


    