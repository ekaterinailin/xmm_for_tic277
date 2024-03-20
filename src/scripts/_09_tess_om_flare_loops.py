import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.constants import mu0, m_e

from astropy.table import Table
import paths

from scipy.optimize import curve_fit


def exponential_decay(t, t0, tau, ampl):

    ex = ampl * np.exp(-(t - t0) / tau)

    ex[t < t0] = 0

    return  ex

if __name__ == "__main__":

    # define sectors
    sectors = 12, 37, 39, 64, 65

    # read in light curves
    lcrs = [Table.read(paths.data / f"tic277539431_tess_detrended_{s}.fits") for s in sectors]

    # read in the flare table
    df = pd.read_csv(paths.data / "tess_flares.csv")

    # delete the largest flare
    df = df.drop(df[df.ed_rec == df.ed_rec.max()].index)

    # select the columns we want, sort by time
    sel = df[["tstart", "tstop", "ampl_rec", 'ed_rec', 'ed_rec_err', 'Sector']].sort_values("tstart")

    efolds = []

    # loop over axes, light curves and sectors
    for lcr, sector in zip(lcrs, sectors):

        # select all flares in this sector
        sel_sector = sel[sel.Sector == sector]

        # mask all flares and calculate the median
        mask = np.ones(len(lcr), dtype=bool)
        for i, flare in sel_sector.iterrows():
            mask &= (lcr['TIME'] < flare.tstart - 0.03) | (lcr['TIME'] > flare.tstop + 0.03)

        # calculate the median
        median = np.median(lcr['DETRENDED_FLUX'][mask])

        # divide by the median
        lcr['DETRENDED_FLUX'] /= median


        # loop over flares
        for i, flare in sel_sector.iterrows():

            # select the light curve
            lc = lcr[(lcr['TIME'] > (flare.tstart - 0.03)) & (lcr['TIME'] < (flare.tstop + 0.03))]

            # select the flare
            lcf = lcr[(lcr['TIME'] > (flare.tstart)) & (lcr['TIME'] < (flare.tstop))]

            # get peak time within tstart and tstop
            peak_time = lcf['TIME'][np.argmax(lcf['DETRENDED_FLUX'])] 


            # rename time and flux
            t, y = lc['TIME'], lc['DETRENDED_FLUX']

            # assert that the x0 values are within the bounds
            assert flare.tstart < peak_time < flare.tstop
            assert 0 < (flare.tstop-flare.tstart)/2 < 1    

            # fit exponential decay with curve_fit
            popt, pcov = curve_fit(exponential_decay, t, y-1,
                                p0=[peak_time, (flare.tstop-flare.tstart)/2, flare.ampl_rec],      
                                bounds=([flare.tstart-0.99652777778 , 0, flare.ampl_rec * 0.95],
                                        [flare.tstop+1.00347222222 , 1, flare.ampl_rec * 1.00347]))
            
            efolds.append(popt[1])

            plt.figure()

            # plot the fit
            plt.plot(t, exponential_decay(t, *popt) + 1, label=f"sector {sector} fit")
            

            # plot the light curve
            plt.plot(lc['TIME'], lc['DETRENDED_FLUX'], label=f"sector {sector}")

            # plot the flare
            plt.axvline(flare.tstart, color='r', ls='--')
            plt.axvline(flare.tstop, color='r', ls='--')

            plt.legend()
            plt.savefig(paths.figures / f"expfit_flare_{sector}_{i}.png")
            plt.close()

    df["efold"] = efolds


    # plot B and L relations

    energies = np.logspace(29.3, 36, 100)
    times = np.logspace(-2, 3, 100) 
    radius = (0.145 * u.R_sun).to(u.cm).value
    E0, t0, B0, L0 = 1.5e30, 3.5, 57, 2.4e9

    for B in [30,60,100, 250, 450]:
        t = t0 * (energies / E0)**(1/3) * (B / B0)**(-5/3) 

        plt.plot(energies , t, linestyle='--', c="grey", alpha=0.8)
        plt.text(energies[-20], t[-22], f"${B}$ G", ha='left', va='center', rotation=20)    

    for L in [5e8,1e9, 3e9, 5e9, 1e10, 2e10, 5e10, 1e11]:
        Lrstar = L / radius
        print(fr"Loop length: {L:.0e} or {Lrstar:.2e} R_*")
        t = t0 * (energies / E0)**(-0.5) * (L / L0)**2.5 
        
        plt.plot(energies , t, linestyle=':', c="grey")
        # covert L to latex
        latexL = f"{L:.0e}".split("e+")[1].lstrip("0")
        latex1 = f"{L:.0e}".split("e+")[0]
        if latex1 == "1":
            latexL = r"10^{" + latexL + r"}"
        else:
            latexL = latex1 + r"\cdot 10^{" + latexL + r"}"

        # find a good position for the annotation
        i = np.where((t < times[-1]))[0][0]
    
        if L > 1e10:
            plt.text(energies[i+8], t[i+10], f"${latexL}$ cm", ha='left', va='center', rotation=-30)
        else:
            plt.text(energies[i], t[i], f"${latexL}$ cm", ha='left', va='center', rotation=-30)

    plt.scatter(df.ed_rec, df.efold * 24 * 60, c="olive", alpha=0.8, s=40)
    plt.scatter([3.7e30],[24/60], c="r", marker="X", label="TIC 277: OM flare", s=40)
    plt.scatter([10**34.473],[0.05*24*60], c="olive", marker='d', s=40, label="TIC 277: TESS high latitude flare")


    handle = plt.scatter([],[], c="olive", label="TIC 277: other TESS flares")

    # add in the Maehara + 2021 flares from YZ CMi
    yzcmi = pd.read_csv(paths.data / "maehara2021_yz_cmi.csv", skiprows=8)
    
    # replace "na" with np.nan
    yzcmi = yzcmi.replace("na", np.nan)
    yzcmi = yzcmi.dropna()
    yzcmi["Ebol(erg)"] = yzcmi["Ebol(erg)"].astype(float)
    yzcmi["eftime(min)"] = yzcmi["eftime(min)"].astype(float)


    plt.scatter(yzcmi["Ebol(erg)"], yzcmi["eftime(min)"], c="grey", marker='.',
                label="M4 YZ CMi: TESS flares (Maehara+2021)", s=40, alpha=0.3, zorder=-10)
    
    # ramsay 2021
    efold = [21.2, 16.1, 25.2, 31.3, 20.2, 11.9, 19.5, 39, 16.6, 19.8, 32.1]
    energy = [6e35, 1e34, 7.5e35, 6.1e35, 8.2e34, 1.1e34, 3.4e35, 4e34, 1.6e34, 4.3e35, 1.3e34]

    ramsay = pd.DataFrame({"Ebol(erg)": energy, "eftime(min)": efold})

    plt.scatter(ramsay["Ebol(erg)"], ramsay["eftime(min)"], c="grey", marker='x',
                label="M3-M4 dwarf: TESS flares (Ramsay+2021)", s=40, alpha=1, zorder=-10)


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("$E_{flare}$ [erg]", fontsize=13)

    plt.ylabel(r"Flare decay $e$-folding time [min]", fontsize=13)
    plt.xlim(energies[0], energies[-1])
    plt.ylim(times[0], times[-1])

    # add legend handle for olive TESS flare circles
    plt.legend(loc=4, frameon=True, fontsize=11)
    plt.tight_layout()
    plt.savefig(paths.figures / "tess_flares_B_L_relation.png")
