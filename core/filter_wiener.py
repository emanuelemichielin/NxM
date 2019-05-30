import numpy as np
from matplotlib import pyplot as plt

def filter_wiener( S, J ):
    S_fft = np.fft.fft(S).tolist()

    for n in range(len(S_fft)):
        if S_fft[ n ]==0: continue
        S_fft[ n ] *= ( 1.0 - ( J[ n ]/( np.absolute( S_fft[ n ] )**2 ) ) )

    return [ coef.real for coef in np.fft.ifft( S_fft ).tolist() ]
