import sys
sys.path.append("../core/")
import math
import numpy as np
import os
from matplotlib import pyplot as plt
from scipy.signal import find_peaks, wiener, peak_widths

from ROOT import TFile, TGraph, TH1F

from DataReader           import DataReader
from NoisePSDGenerator    import NoisePSDGenerator
from save_noise_psds      import save_noise_psds
from get_noise_psds       import get_noise_psds
from calc_r               import calc_r
from calc_theta           import calc_theta
from estimate_bins_part   import estimate_bins_part
from TemplateGeneratorNxM import TemplateGeneratorNxM
from save_templates_nxm   import save_templates_nxm
from get_templates_nxm    import get_templates_nxm
from OFManagerNxM         import OFManagerNxM
from filter_wiener        import filter_wiener

integral = []

# path to directory containing data files
filepath='/gpfs/slac/staas/fs1/g/supercdms//data/CDMS/SLAC/R56/Raw/09190602_1927'
 
# specifies series to be analyzed
series=["09190602_1927_F00"+str(i)+".mid.gz" for i in range(70,71)] 


dr = DataReader()
dr.OpenFile( filepath, series, 0 )

V = get_noise_psds( 'noise_psds.gz' )

while dr.LoadEvent(trigger='Trigger'):
    S = dr.GetTraces()
    noise_w = np.array([V[a][a] for a in range(len(S))])
    dataS = np.asarray(S)
    dataS = np.sum(dataS, axis=0)
    if (np.mean(dataS[0:3000])>np.mean(dataS[3000:10000]) and np.mean(dataS[0:3000])>np.mean(dataS[10000:15000]) and np.mean(dataS[0:3000])>np.mean(dataS[30000:32000])): continue
    if (np.mean(dataS[0:3000]) > 1.05*np.mean(dataS[30000:32000])) : continue 
    if (1.05*np.mean(dataS[0:3000]) < np.mean(dataS[30000:32000])) : continue
    dataS = dataS*1E7
    peaks, properties = find_peaks(dataS.transpose(), prominence=1, width=20)
    width_half = peak_widths(dataS.transpose(), peaks, rel_height=0.5)
    
    if (len(peaks) == 0 or len(peaks) > 1 or any(width_half)>2000) :
        continue
    tmp_S = np.array(S)
    avg =  np.mean(tmp_S[:,0:5000],axis=1) 
    for i in range(len(S)):
        S[i] = S[i]-avg[i]
        
    amps = [ sum( wiener( S[ a ], mysize=75, noise = noise_w[a].real ) ) for a in range(len(S)) ]
    integral.append(sum( amps ))

plt.hist(integral, 200, range=(0,0.1))
plt.show()
