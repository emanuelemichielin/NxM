import math
import numpy as np

from utils import interpolate_parab

class OptimalFilter:

    def __init__( self, dt, t_pre, U, J ):

        # dt is width of the time bins
        # t_pre is the interval of time between the start of the trace and the trigger
        # Note that t_pre is positive
        # U is the template
        # J is the PSD of the noise

        self.dt = dt
        self.t_pre = t_pre
        self.n_trig = int( t_pre/dt )

        self.num_bins_t = len( U )
        self.loop_bins_t = xrange( self.num_bins_t )

        self.U = [ U[ n ] for n in self.loop_bins_t ]
        U_fft =  np.fft.fft( U ).tolist()

        self.J = [ J[ n ] for n in self.loop_bins_t ]

        self.denom = 0.0

        for n in self.loop_bins_t:
            self.denom += ( np.absolute( U_fft[ n ] )**2 )/J[ n ].real

        self.F = []

        for n in self.loop_bins_t:
            self.F.append( np.conj( U_fft[ n ] )/( self.denom*J[ n ] ) )

        self.sigma = 1.0/math.sqrt( self.denom )

        self.S_fft = None

        self.result = { 'A0'    : None,
                        'chisq0': None,
                        't0'    : None,
                        'A'     : None,
                        'chisq' : None }

    def GetTheoreticalResolution( self ):
        return self.sigma

    def Execute( self, S ):
        self.S_fft = np.fft.fft( S ).tolist()

        func_a = np.real( np.fft.ifft( [ self.F[ n ]*self.S_fft[ n ] for n in self.loop_bins_t ] ) ).tolist()

        for n in self.loop_bins_t:
            func_a[ n ] *= float( self.num_bins_t )

        self.result[ 'A0' ] = float( func_a[ 0 ] )

        chisq_base = 0.0

        for n in self.loop_bins_t:
            chisq_base += ( np.absolute( self.S_fft[ n ] )**2 )/self.J[ n ]

        chisq_base = chisq_base.real

        self.result[ 'chisq0' ] = float( chisq_base-( func_a[ 0 ]**2 )/self.denom )

        n_min = func_a.index( max( func_a ) )

        if n_min >= self.num_bins_t-self.n_trig:
            n_min -= self.num_bins_t

        if n_min > -self.n_trig and n_min < ( self.num_bins_t-self.n_trig )-1:
            n_prev = n_min-1
            n_post = ( n_min+1 )%self.num_bins_t

            parab_a = interpolate_parab( func_a[ n_prev ], func_a[ n_min ], func_a[ n_post ] )
            dn0 = -parab_a[ 1 ]/( 2.0*parab_a[ 0 ] )

            self.result[ 't0' ] = float( ( 0.5+float( n_min )+dn0 )*self.dt )
            self.result[ 'A'  ] = float( ( parab_a[ 0 ]*( dn0**2 ) )+( parab_a[ 1 ]*dn0 )+parab_a[ 2 ] )

            chisq_min = chisq_base-( func_a[ n_min ]**2 )/self.denom
            chisq_prev = chisq_base-( func_a[ n_prev ]**2 )/self.denom
            chisq_post = chisq_base-( func_a[ n_post ]**2 )/self.denom

            parab_chisq = interpolate_parab( chisq_prev, chisq_min, chisq_post )
            self.result[ 'chisq' ] = float( parab_chisq[ 2 ]-( ( parab_chisq[ 1 ]**2 )/( 4.0*parab_chisq[ 0 ] ) ) )

        else:
            self.result[ 't0'    ] = float( ( 0.5+float( n_min ) )*self.dt )
            self.result[ 'A'     ] = float( func_a[ n_min ] )
            self.result[ 'chisq' ] = float( chisq_base-( func_a[ n_min ]**2 )/self.denom )

        return self.result

    def reset_result( self ):
        self.result[ 'A0'     ] = None
        self.result[ 'chisq0' ] = None
        self.result[ 't0'     ] = None
        self.result[ 'A'      ] = None
        self.result[ 'chisq'  ] = None

        return True
