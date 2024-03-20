import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import paths
from scipy.optimize import curve_fit



def exponential_decay(t, t0, tau, ampl):

    ex = ampl * np.exp(-(t - t0) / tau)

    ex[t < t0] = 0

    return  ex


if __name__ == '__main__':

    df = pd.read_csv(paths.data / 'timeseries.csv')

    # calculate the median of the time series
    med = df.median()

    # mask the data to only include the flare and its surroundings
    mask = (df.time > 7.76075e8) & (df.time < 7.76077e8)
    df = df[mask]

    # normalize the data to the median
    t, y = df.time, df.rate / med.rate

    # find the peak and the peak amplitude
    peak = df.time[df.rate==df.rate.max()]
    peaka = df.rate.max() / med.rate - 1

    # fit exponential decay with curve_fit
    popt, pcov = curve_fit(exponential_decay, t, y-1,
                        p0=[peak, 20, peaka],      
                        bounds=([peak-10 , 1, peaka * 0.95],
                                [peak+10 , 100, peaka * 1.00347]))
    
    # plot the data and the fit
    plt.figure()
    plt.plot(t, exponential_decay(t, *popt) + 1)
    plt.plot(t, y)
    plt.xlim(7.760755e8, 7.760765e8)
    plt.xlabel('Time [s]')
    plt.ylabel('Normalized flux')
    plt.savefig(paths.figures / 'OM_exponential_decay.png')

    # print e-folding time
    print('e-folding time [s]:', popt[1])