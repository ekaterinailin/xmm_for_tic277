"""
Python 3.8 - UTF-8

X-ray Loops
Ekaterina Ilin, 2023
MIT License

---

This script reads in high level results and produces latex values in the output
folder.

Values covered:
- FFD slope
- FFD intercept
- R31.5

- Lx
- Fx
- Lx/Lbol

- T1 and T2
- quiescient mean T
- flaring mean T

- flare energy OM and EPIC
"""

import pandas as pd
import numpy as np
import paths

if __name__ == "__main__":
    

    # FFD alpha beta -----------------------------------------------------------

    df = pd.read_csv(paths.data / "tess_ffd.csv")
    row = df[df.TIC == 277539431].iloc[0]

    # convert alpha and beta results to latex strings with 
    # upper and lower uncertainties for beta convert to log10(beta)
    alpha, beta = row["alpha"], row["beta"]
    beta_low_err, beta_up_err = row["beta_low_err"], row["beta_up_err"]

    alpha_str = f"${alpha:.2f}_{{-{row['alpha_low_err']:.2f}}}^{{+{row['alpha_up_err']:.2f}}}$"
    beta_str = f"${np.log10(beta):.2f}_{{-{np.log10(beta_low_err):.2f}}}^{{+{np.log10(beta_up_err):.2f}}}" + r"\,\mathrm{d}^{-1}$"

    with open(paths.output / "tess_ffd_alpha.tex", "w") as f:
        f.write(f"{alpha_str}")

    with open(paths.output / "tess_ffd_beta.tex", "w") as f:
        f.write(f"{beta_str}")

    # R31.5 --------------------------------------------------------------------

    r315 = beta / (alpha - 1) * (10**31.5)**(-alpha +1)
    r315 = np.log10(r315)

    r315str = f"${r315:.2f}" + r"\,\mathrm{d}^{-1}$"

    print(r315str)

    with open(paths.output / "R315.tex", "w") as f:
        f.write(r315str)

    # Lx -----------------------------------------------------------------------

    df = pd.read_csv(paths.data / "joint_vapec_chain_fits.csv")

    # select full data set
    row = df[df.subset == "noflare"].iloc[0]

    # get Lx and make latex string with a single uncertainty
    Lx, Lxerr = row.Lx_erg_s, row.Lx_erg_s_err
    Lx_str = fr"${Lx/1e26:.2f} \pm {Lxerr/1e26:.2f}" + r"\times 10^{26}\,\rm{erg s}^{-1}$"

    print(Lx_str)

    with open(paths.output / "epic_Lx.tex", "w") as f:
        f.write(f"{Lx_str}")

    # Fx -----------------------------------------------------------------------

    flux, fluxerr = row.flux_erg_s_cm2, row.flux_erg_s_cm2_err
    flux_str = fr"${flux*1e14:.1f} \pm {fluxerr*1e14:.1f}" + r"\times 10^{-14} \,\rm{erg}\, \rm{cm}^{-2} \rm{s}^{-1}$"

    print(flux_str)

    with open(paths.output / "epic_flux.tex", "w") as f:
        f.write(f"{flux_str}")

    # L_bol and L_X / L_bol --------------------------------------------------------------

    ilin2021 = pd.read_csv(paths.data / "ilin2021updated_w_Rossby_Lbol.csv")
    row = ilin2021[ilin2021.TIC == 277539431].iloc[0]
    Lbol, eLbol = row.Lbol_erg_s, row.eLbol_erg_s

    Lbolstr = fr"${Lbol/1e30:.1f} \pm {eLbol/1e30:.1f}" + r" \times 10^{30}\,\rm{erg s}^{-1}$"

    with open(paths.output / "Lbol.tex", "w") as f:
        f.write(f"{Lbolstr}")

    lxlbol, elxlbol = (Lx / Lbol, 
                    Lx / Lbol * np.sqrt(Lxerr**2 / (Lx**2) + eLbol**2 / (Lbol**2)))


    lxlbolstr = fr"${lxlbol*1e4:.1f}\pm {elxlbol*1e4:.1f}" + r"\times 10^{-4}$ "

    print(lxlbolstr)

    with open(paths.output / "lxlbol.tex", "w") as f:
        f.write(rf"{lxlbolstr}")

    # T1 and T2 ----------------------------------------------------------------

    df = pd.read_csv(paths.data / "mcmc_results.csv")
    df = df.rename(columns={"Unnamed: 0": "subset"})
    row = df[df.subset == "full data set"].iloc[0]

    # get T1 and make latex string with a upper and lower uncertainty
    # derived from the 16th and 84th percentile of the posterior distribution
    T150, T116, T184 = row.T1_50, row.T1_16, row.T1_84
    T1_str = f"${T150:.1f}_{{-{T150-T116:.1f}}}^{{+{T184-T150:.1f}}}\,$MK"

    with open(paths.output / "T1.tex", "w") as f:
        f.write(f"{T1_str}")

    # get T2 and make latex string with a upper and lower uncertainty
    # derived from the 16th and 84th percentile of the posterior distribution
    T250, T216, T284 = row.T2_50, row.T2_16, row.T2_84
    T2_str = f"${T250:.1f}_{{-{T250-T216:.1f}}}^{{+{T284-T250:.1f}}}\,$MK"

    with open(paths.output / "T2.tex", "w") as f:
        f.write(f"{T2_str}")

    # quiescent mean T ---------------------------------------------------------

    row = df[df.subset == "quiescent"].iloc[0]
    T1q, T2q = row.T1_50, row.T2_50
    norm1q, norm2q = row.norm1_50, row.norm2_50
    Tqmean = (T1q * norm1q + T2q * norm2q) / (norm1q + norm2q)
    Tqmean_str = f"${Tqmean:.1f}\,$MK"

    with open(paths.output / "Tqmean.tex", "w") as f:
        f.write(f"{Tqmean_str}")

    # flaring mean T------------------------------------------------------------

    row = df[df.subset == "flaring"].iloc[0]
    T1f, T2f = row.T1_50, row.T2_50
    norm1f, norm2f = row.norm1_50, row.norm2_50
    Tfmean = (T1f * norm1f + T2f * norm2f) / (norm1f + norm2f)
    Tfmean_str = f"${Tfmean:.1f}\,$MK"

    with open(paths.output / "Tfmean.tex", "w") as f:
        f.write(f"{Tfmean_str}")

    # flare energy OM and EPIC -------------------------------------------------

    epicom = pd.read_csv(paths.data / "flare_energies.csv")
    epic = epicom[epicom.instrument == "EPIC"].iloc[0]
    om = epicom[epicom.instrument == "OM"].iloc[0]
    E_epic, eE_epic = epic.E_erg, epic.eE_erg
    E_om, eE_om = om.E_erg, om.eE_erg

    e, ee = E_epic / 1e30, eE_epic /1e30
    epicstr = fr"${e:.1f}\pm{ee:.1f}" + r"\times 10^{30}\,$erg"

    with open(paths.output / "epic_flare.tex", "w") as f:
        f.write(epicstr)

    e, ee = E_om / 1e30, eE_om /1e30
    omstr = fr"${e:.1f}\pm{ee:.1f}" + r"\times 10^{30}\,$erg"

    with open(paths.output / "om_flare.tex", "w") as f:
        f.write(omstr)

