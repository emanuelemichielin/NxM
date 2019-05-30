from scdmsPyTools.BatTools.IO import *
import sys
sys.path.append("../core/")
from pandas import * 
import numpy as np 
from pylab import *
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks, peak_widths, wiener
from calc_r               import calc_r
from calc_theta           import calc_theta

# path to directory containing data files
filepath='/nfs/slac/g/supercdms/data/CDMS/SLAC/R52/Raw/09190411_2241'
 
# specifies series to be analyzed
#series="09190411_2241"
series=["09190411_2241_F0030.mid.gz"]  

# loads raw data
ev = getRawEvents(filepath,series, outputFormat=2)


channels1 = ['PAS1','PBS1','PCS1','PDS1','PES1','PFS1']
channels2 = ['PAS2','PBS2','PCS2','PDS2','PES2','PFS2']

 
# specify detector and channel
for k in range(0,20):
    dataS1=[]
    dataS2=[]

    for i in channels1:
        dataS1.append(ev[k]['Z6'][i])
    for i in channels2:
        dataS2.append(ev[k]['Z6'][i])




    dataS1 = np.asarray(dataS1)
    dataS2 = np.asarray(dataS2)
    data_sumS1 = np.sum(dataS1, axis=0)
    data_sumS2 = np.sum(dataS2, axis=0)

    #if (np.mean(data_sumS1[0:3000])>np.mean(data_sumS1[3000:10000]) and np.mean(data_sumS1[0:3000])>np.mean(data_sumS1[10000:15000]) and np.mean(data_sumS1[0:3000])>np.mean(data_sumS1[30000:32000])): continue
    #if (np.mean(data_sumS1[0:3000]) > 1.05*np.mean(data_sumS1[30000:32000])) : continue 
    #if (1.05*np.mean(data_sumS1[0:3000]) < np.mean(data_sumS1[30000:32000])) : continue

    peaks, properties = find_peaks(data_sumS1.transpose(), prominence=1, width=20)
    width_half = peak_widths(data_sumS1.transpose(), peaks, rel_height=0.5)

    if (len(peaks) > 1 or any(width_half)>4000) :
        continue
    tmp_S = np.array(dataS1)
    avg =  np.mean(tmp_S[:,0:5000],axis=1) 
    for i in range(len(dataS1)):
        dataS1[i] = dataS1[i]-avg[i]

    tmp_S = np.array(dataS2)
    avg =  np.mean(tmp_S[:,0:5000],axis=1) 
    for i in range(len(dataS2)):
        dataS2[i] = dataS2[i]-avg[i]



    amps_S1 = [ sum( wiener( dataS1[ a ], mysize=75) ) for a in range(len(dataS1)) ]
    amps_S2 = [ sum( wiener( dataS2[ a ], mysize=75) ) for a in range(len(dataS2)) ]
        
    theta_S1 =  calc_theta( amps_S1 ) 
    r_S1 = calc_r( amps_S1 ) 

    theta_S2 =  calc_theta( amps_S2 ) 
    r_S2 = calc_r( amps_S2 ) 


    n = 32768
    r = 625000
    x = np.linspace(0, 0.0524288, 32768)
    freq = np.fft.fftfreq(n, d = 1/r)
    k = freq > 0

    plt.figure(0, dpi=120)  
    for i,chan in enumerate(channels1):

        sp1 = np.fft.fft(dataS1[i].transpose())

        plt.subplot(111)
        plt.xlabel("Time [s]")
        plt.plot(x, dataS1[i].transpose(), label= chan, alpha=0.8)
        plt.xlim(0.024,0.035)
        plt.subplot(212)
        plt.loglog(freq[k], np.abs(sp1[k]) ** 2, label= chan, alpha=0.8)
        plt.xlabel("Frequency [Hz]")
        #if (i==5):
        #    peaks, properties = find_peaks( np.log(np.abs(sp1)**2), prominence=10)
        #    print(freq[peaks])
        #    plt.loglog((freq[peaks]), ((np.abs(sp1[peaks])**2)), "x")
    plt.grid()
    plt.legend()
    plt.savefig("png/nice_zoomS1.png")



    plt.figure(1, dpi=120)  
    for i,chan in enumerate(channels2):

        sp2 = np.fft.fft(dataS2[i].transpose())

        plt.subplot(111)
        plt.xlabel("Time [s]")
        plt.plot(x,dataS2[i].transpose(), label= chan, alpha=0.8)
        plt.xlim(0.024,0.035)
        #plt.subplot(212)
        #plt.xlabel("Frequency [Hz]")   
        #plt.loglog(freq[k], np.abs(sp2[k]) ** 2, label= chan, alpha=0.8 )
 
    plt.legend()
    plt.grid()
    plt.savefig("png/nice_zoomS2.png")



    plt.show()
