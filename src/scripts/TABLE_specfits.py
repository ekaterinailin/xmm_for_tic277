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

    # read in the data
    df = pd.read_csv(paths.data / "apec_apec_results.csv")
    df["model"] = "APEC+APEC"
    print(df.shape)
    dfjoint0 = pd.read_csv(paths.data / "vapec_vapec_results.csv")
    dfjoint0["model"] = "VAPEC+VAPEC"
    dfjoint = pd.read_csv(paths.data / "apec_apec_joint.csv")
    dfjoint["model"] = "APEC+APEC"
    dfjoint2 = pd.read_csv(paths.data / "vapec_vapec_joint.csv")
    dfjoint2["model"] = "VAPEC+VAPEC"

    print(dfjoint.columns)
    # append the joint results
    df = pd.concat([df, dfjoint0, dfjoint, dfjoint2], ignore_index=True)
    print(df.columns)

    # values with one error
    cols2 = [("Lx_erg_s","Lx_erg_s_err"),
            ("Lx_Lbol","e_Lx_Lbol") ]

    # values with up and lower 90% errors
    cols_mid_low_high = [("flux_erg_s_cm2","flux_erg_s_cm2_low","flux_erg_s_cm2_high"),
                        ("T_MK_1","T_MK_1_low","T_MK_1_high"),
                        ("T_MK_5","T_MK_5_low","T_MK_5_high"),]
    
    # rename columns
    mapcolnames = {"flux_erg_s_cm2":r"$10^{-14} F_X$ [erg/s/cm$^2$]",
                    "T_MK_1":r"$T_1$ [MK]",
                    "T_MK_5":r"$T_2$ [MK]",
                    "Lx_erg_s":r"$10^{26} L_X$ [erg/s]",
                    "Lx_Lbol":r"$10^{-4} L_X/L_{bol} $"}

    # convert to new units
    df["Lx_Lbol"] = 1e4 * df["Lx_Lbol"]
    df["e_Lx_Lbol"] = 1e4 * df["e_Lx_Lbol"]
    df["Lx_erg_s"] = 1e-26 * df["Lx_erg_s"]
    df["Lx_erg_s_err"] = 1e-26 * df["Lx_erg_s_err"]
    df["flux_erg_s_cm2"] = 1e14 * df["flux_erg_s_cm2"]
    df["flux_erg_s_cm2_low"] = 1e14 * df["flux_erg_s_cm2_low"]
    df["flux_erg_s_cm2_high"] = 1e14 * df["flux_erg_s_cm2_high"]

    # convert to latex strings
    for col, err in cols2:
     
        newname = mapcolnames[col]
        df[newname] = tex_one_err(df[col], df[err])
        del df[col]
        del df[err]

    # convert to latex strings
    for col, low, up in cols_mid_low_high:
        
        newname = mapcolnames[col]
        df[newname] = tex_up_low(df[col], df[up]-df[col], df[low]-df[col])
        del df[col]
        del df[up]
        del df[low]


    # delete xtra columns
    del df["detector"]
    del df["cut"]

    for col in df.columns:
        if "Rossby" in col:
            del df[col]


    # convert to LaTeX
    df = df.sort_values("model")
    string = df.to_latex(escape=False,index=False)

    # layout
    string = string.replace("midrule","hline")
    string = string.replace("toprule","hline")
    string = string.replace("bottomrule","hline")
    string = string.replace("label","Detector (Cut)")
    string = string.replace("$nan[nan]$","-")
    string = string.replace(r"$nan^{nan}_{nan}$","-")

    # write to file
    path = paths.output / f"specfit.tex"
    print("Write to file: ", path)

    with open(path, "w") as f:
        f.write(string)