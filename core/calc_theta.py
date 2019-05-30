import math
import numpy as np

def calc_theta( A ):
    x_u_theta = [ 1.0, -0.5                , -0.5                  ]
    y_u_theta = [ 0.0, 0.5*math.sqrt( 3.0 ), -0.5*math.sqrt( 3.0 ) ]

    return math.atan2( np.average( y_u_theta, weights = A[ 2:5 ] ), np.average( x_u_theta, weights = A[ 2:5 ] ) )
