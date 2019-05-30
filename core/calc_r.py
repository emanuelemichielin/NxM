import numpy as np

def calc_r( A ):
    return np.average( [ 0.0, 4.37 ], weights = [ A[ 5 ], 0.5*sum( A[ 0:2 ] ) ] )
