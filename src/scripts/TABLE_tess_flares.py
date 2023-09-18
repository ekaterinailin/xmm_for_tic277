"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in the TESS flare table and produces a LaTeX table.
"""

import pandas as pd
import paths

if __name__ == "__main__":

    # read in the flare table
    df = pd.read_csv(paths.data / "tess_flares.csv")

    # select the columns we want, sort by time
    sel = df[["tstart", "ampl_rec", 'ed_rec', 'ed_rec_err', 'Sector']].sort_values("tstart")

    # convert to 10^30 erg
    sel["ed_rec"] = sel["ed_rec"] / 1e31
    sel["ed_rec_err"] = sel["ed_rec_err"] / 1e31

    # format
    sel[r"$E_{\rm bol}$ [$10^{31}$ erg]"] = (sel["ed_rec"].round(1).astype(str) + 
                                            " [" + 
                                            sel["ed_rec_err"].round(1).astype(str) + 
                                            "]")
    # delete unused columns
    del sel["ed_rec_err"]
    del sel["ed_rec"]

    # rename columns
    sel = sel.rename(columns={"tstart": r"$t_{s}$ [BJD - 2457000]", "ampl_rec": r"$a$"})

    # select columns in the correct order
    sel = sel[[r"$t_{s}$ [BJD - 2457000]", r"$a$", r"$E_{\rm bol}$ [$10^{31}$ erg]", "Sector"]]

    # convert to LaTeX
    string = sel.to_latex(index=False, escape=False)

    # layout
    string = string.replace("midrule","hline")
    string = string.replace("toprule","hline")
    string = string.replace("bottomrule","hline")

    # write to file
    with open(paths.output / "tess_flares.tex", "w") as f:
        f.write(string)