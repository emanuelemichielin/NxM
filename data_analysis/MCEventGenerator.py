import math
import numpy as np

from ROOT import TGraph, TH1F, TRandom3

from utils import create_name_hist, get_angle_std, check_angle, vector_distribution, draw_hists, draw_graphs

class MCEventGenerator:

    def __init__( self, seed, correlated = False ):
        num_bins_t = 1024
        dt = 0.0016
        t_pre = 0.2655

        num_bins_margin_t = 50

        num_bins_l = 50
        l_max = 10.0

        self.max_radius = 3.5
        self.sigma_noise = 0.0002
        self.threshold_trig = 0.0005

        self.list_r_bot = [ 0.0   , 1.9596, 1.9596, 1.9596, 3.9192, 4.3818 ]
        self.list_r_top = [ 1.9596, 3.9192, 3.9192, 3.9192, 4.3818, 4.8    ]
        self.list_theta_clockw = [ -math.pi, -math.pi/3.0, math.pi/3.0, -math.pi    , -math.pi, -math.pi ]
        self.list_theta_anticw = [ math.pi , math.pi/3.0 , math.pi    , -math.pi/3.0, math.pi , math.pi  ]

        self.initialized = False

        self.n_trig = int( t_pre/dt )

        if self.n_trig < num_bins_t and self.n_trig > num_bins_margin_t and num_bins_t-self.n_trig > num_bins_margin_t:
            self.num_channels = len( self.list_r_bot )
            self.num_bins_t = num_bins_t

            self.num_bins_margin_t = num_bins_margin_t

            self.loop_channels = xrange( self.num_channels )
            self.loop_bins_t = xrange( self.num_bins_t )

            self.list_t_n = []

            for n in xrange( self.n_trig, num_bins_t+num_bins_margin_t ):
                self.list_t_n.append( ( ( 0.5+float( n ) )*dt )-t_pre )

            dl = l_max/float( num_bins_l )
            self.loop_bins_l = xrange( num_bins_l )

            self.list_l_i = [ float( i )*dl for i in self.loop_bins_l ]
            self.list_c_i = [ [ None for i in self.loop_bins_l ] for a in self.loop_channels ]
            self.list_f_i = [ [ self.calc_model( t_n, l_i ) for t_n in self.list_t_n ] for l_i in self.list_l_i ]

            self.update_list_c_i( 1.0, 0.0, 0.0 )

            norm = 0.0

            for a in self.loop_channels:
                for i in self.loop_bins_l:
                    for n in xrange( len( self.list_t_n ) ):
                        norm += self.list_c_i[ a ][ i ]*self.list_f_i[ i ][ n ]

            for i in self.loop_bins_l:
                for n in xrange( len( self.list_t_n ) ):
                    self.list_f_i[ i ][ n ] /= norm

            self.rand = TRandom3( seed )

            if correlated:
                self.phases_ref = None
            
            self.correlated = correlated

            self.S = None
            self.last_mc_truth = None

            self.vector_distributions = [ vector_distribution() for a in self.loop_channels ]

            self.initialized = True

    def calc_model( self, t, l ):
        if t < 0.0:
            return 0.0

        a_T = 0.01
        b_T = 0.05
        eta0 = 2.0
        l0 = 1.0
        tau1 = 0.1
        tau2 = 0.2

        T = a_T+( b_T*l )
        eta = eta0*( ( 1.0+( T/tau1 ) )/( 1.0+( T/tau2 ) ) )*math.exp( -l/l0 )

        return ( 1.0-math.exp( -t/T ) )*( ( eta*math.exp( -t/tau1 ) )+( ( l/l0 )*math.exp( -t/tau2 ) ) )

    def GetNumChannels( self ):
        return self.num_channels

    def GetNumBinsT( self ):
        return self.num_bins_t

    def GenerateNoiseTraces( self ):
        if not self.initialized:
            return None

        if self.correlated:
            self.phases_ref = [ self.rand.Uniform( -math.pi, math.pi ) for n in xrange( ( self.num_bins_t+2 )/2 ) ]

        S = [ self.generate_noise_trace() for a in self.loop_channels ]

        if type( self.last_mc_truth ) == tuple:
            self.last_mc_truth = None

        return S

    def generate_noise_trace( self ):
        coefs_fft = []

        for n in xrange( ( self.num_bins_t+2 )/2 ):
            if self.correlated:
                phase = self.phases_ref[ n ]+math.asin( self.rand.Uniform( -1.0, 1.0 ) )
                coefs_fft.append( math.cos( phase )+( 1j*math.sin( phase ) ) )

            else:
                phase = self.rand.Uniform( -math.pi, math.pi )
                coefs_fft.append( math.cos( phase )+( 1j*math.sin( phase ) ) )

        return [ self.sigma_noise*coef for coef in np.fft.irfft( coefs_fft, self.num_bins_t ).tolist() ]

    def GenerateTraces( self, list_par = [] ):
        if not self.initialized:
            return None

        if len( list_par ) == 0:
            E_gen = 1.0
            r_gen = self.max_radius*math.sqrt( self.rand.Uniform() )
            theta_gen = ( 2.0*math.pi*self.rand.Uniform() )-math.pi

            self.update_list_c_i( E_gen, r_gen, theta_gen )

        elif len( list_par ) == 3:
            E_gen = list_par[ 0 ]
            r_gen = list_par[ 1 ]
            theta_gen = get_angle_std( list_par[ 2 ] )

            self.update_list_c_i( E_gen, r_gen, theta_gen )

        else:
            return None

        pulses = []

        for a in self.loop_channels:
            pulses.append( [] )

            for n in xrange( len( self.list_t_n ) ):
                pulses[ -1 ].append( 0.0 )

                for i in self.loop_bins_l:
                    pulses[ -1 ][ -1 ] += self.list_c_i[ a ][ i ]*self.list_f_i[ i ][ n ]

        if self.correlated:
            self.phases_ref = [ self.rand.Uniform( -math.pi, math.pi ) for n in xrange( ( self.num_bins_t+2 )/2 ) ]

        self.S = []

        shift_n = self.num_bins_margin_t
        status_trig = []

        for a in self.loop_channels:
            self.S.append( self.generate_noise_trace() )

            for n in xrange( self.num_bins_t-self.n_trig ):
                self.S[ -1 ][ self.n_trig+n ] += pulses[ a ][ n ]

            status_trig.append( True )

            for n in xrange( self.num_bins_margin_t ):
                if self.S[ -1 ][ self.n_trig+n ] > self.threshold_trig:
                    if n < shift_n:
                        shift_n = n

                    break

            else:
                status_trig[ -1 ] = False

        if not True in status_trig:
            return None

        if shift_n > 0:
            for a in self.loop_channels:
                for n in xrange( shift_n ):
                    self.S[ a ].append( self.S[ a ].pop( 0 )+pulses[ a ][ ( self.num_bins_t-self.n_trig )+n ] )

        A = [ sum( pulse ) for pulse in pulses ]
        self.last_mc_truth = ( [ E_gen, r_gen, theta_gen ], A, shift_n, [] )

        for a in self.loop_channels:
            self.last_mc_truth[ 3 ].append( [] )

            for n in xrange( self.n_trig-shift_n ):
                self.last_mc_truth[ 3 ][ -1 ].append( 0.0 )

            for n in xrange( ( self.num_bins_t-self.n_trig )+shift_n ):
                self.last_mc_truth[ 3 ][ -1 ].append( pulses[ a ][ n ] )

        return self.S

    def update_list_c_i( self, E_gen, r_gen, theta_gen_, draw = False ):
        theta_gen = get_angle_std( theta_gen_ )

        x0 = r_gen*math.cos( theta_gen )
        y0 = r_gen*math.sin( theta_gen )

        for a in self.loop_channels:
            for i in self.loop_bins_l:
                l = self.list_l_i[ i ]

                list_phi_cut = []

                if i > 0:
                    if r_gen > 0.0:
                        if self.list_r_bot[ a ] > 0.0:
                            arg = ( ( self.list_r_bot[ a ]*self.list_r_bot[ a ] )-( r_gen*r_gen )-( l*l ) )/( 2.0*r_gen*l )
        
                            if math.fabs( arg ) < 1.0:
                                for phi_cut in [ get_angle_std( theta_gen-math.acos( arg ) ),
                                                 get_angle_std( theta_gen+math.acos( arg ) ) ]:
                                    x_cut = x0+( l*math.cos( phi_cut ) )
                                    y_cut = y0+( l*math.sin( phi_cut ) )
                                    theta_cut = math.atan2( y_cut, x_cut )

                                    if check_angle( theta_cut, self.list_theta_clockw[ a ], self.list_theta_anticw[ a ] ):
                                        list_phi_cut.append( phi_cut )

                        arg = ( ( self.list_r_top[ a ]*self.list_r_top[ a ] )-( r_gen*r_gen )-( l*l ) )/( 2.0*r_gen*l )
        
                        if math.fabs( arg ) < 1.0:
                            for phi_cut in [ get_angle_std( theta_gen-math.acos( arg ) ),
                                             get_angle_std( theta_gen+math.acos( arg ) ) ]:
                                x_cut = x0+( l*math.cos( phi_cut ) )
                                y_cut = y0+( l*math.sin( phi_cut ) )
                                theta_cut = math.atan2( y_cut, x_cut )

                                if check_angle( theta_cut, self.list_theta_clockw[ a ], self.list_theta_anticw[ a ] ):
                                    list_phi_cut.append( phi_cut )

                    if self.list_theta_clockw[ a ] != -math.pi or self.list_theta_anticw[ a ] != math.pi:
                        arg = math.sin( theta_gen-self.list_theta_clockw[ a ] )*r_gen/l

                        if math.fabs( arg ) < 1.0:
                            for phi_cut in [ get_angle_std( self.list_theta_clockw[ a ]-math.asin( arg )         ),
                                             get_angle_std( self.list_theta_clockw[ a ]+math.asin( arg )-math.pi ) ]:
                                x_cut = x0+( l*math.cos( phi_cut ) )
                                y_cut = y0+( l*math.sin( phi_cut ) )
                                r_cut = math.sqrt( ( x_cut*x_cut )+( y_cut*y_cut ) )
                                theta_cut = math.atan2( y_cut, x_cut )

                                if math.fabs( math.fabs( theta_cut-self.list_theta_clockw[ a ] )-math.pi ) < 1.0e-2:
                                    continue

                                if r_cut >= self.list_r_bot[ a ] and r_cut < self.list_r_top[ a ]:
                                    list_phi_cut.append( phi_cut )

                        arg = math.sin( theta_gen-self.list_theta_anticw[ a ] )*r_gen/l

                        if math.fabs( arg ) < 1.0:
                            for phi_cut in [ get_angle_std( self.list_theta_anticw[ a ]-math.asin( arg )         ),
                                             get_angle_std( self.list_theta_anticw[ a ]+math.asin( arg )-math.pi ) ]:
                                x_cut = x0+( l*math.cos( phi_cut ) )
                                y_cut = y0+( l*math.sin( phi_cut ) )
                                r_cut = math.sqrt( ( x_cut*x_cut )+( y_cut*y_cut ) )
                                theta_cut = math.atan2( y_cut, x_cut )

                                if math.fabs( math.fabs( theta_cut-self.list_theta_anticw[ a ] )-math.pi ) < 1.0e-2:
                                    continue

                                if r_cut >= self.list_r_bot[ a ] and r_cut < self.list_r_top[ a ]:
                                    list_phi_cut.append( phi_cut )

                arc = 0.0

                if len( list_phi_cut ) == 0:
                    x_test = x0+l
                    y_test = y0
                    r_test = math.sqrt( ( x_test*x_test )+( y_test*y_test ) )
                    theta_test = math.atan2( y_test, x_test )

                    if r_test >= self.list_r_bot[ a ] and r_test < self.list_r_top[ a ]:
                        if check_angle( theta_test, self.list_theta_clockw[ a ], self.list_theta_anticw[ a ] ):
                            arc += 2.0*math.pi

                else:
                    list_phi_cut.sort()

                    phi_cut_prev = list_phi_cut[ -1 ]-( 2.0*math.pi )
                    inside = False

                    phi_test = get_angle_std( 0.5*( list_phi_cut[ 0 ]+phi_cut_prev ) )
                    x_test = x0+( l*math.cos( phi_test ) )
                    y_test = y0+( l*math.sin( phi_test ) )
                    r_test = math.sqrt( ( x_test*x_test )+( y_test*y_test ) )
                    theta_test = math.atan2( y_test, x_test )

                    if r_test >= self.list_r_bot[ a ] and r_test < self.list_r_top[ a ]:
                        if check_angle( theta_test, self.list_theta_clockw[ a ], self.list_theta_anticw[ a ] ):
                            inside = True

                    for phi_cut in list_phi_cut:
                        if inside:
                            arc += phi_cut-phi_cut_prev
                            inside = False

                            if draw:
                                self.update_vector_distributions( a, x0, y0, l, phi_cut_prev, phi_cut )

                        else:
                            inside = True

                        phi_cut_prev = phi_cut

                self.list_c_i[ a ][ i ] = E_gen*arc/( 2.0*math.pi )

    def update_vector_distributions( self, a, x0, y0, l, phi_cut_prev, phi_cut ):
        max_step_size = 0.1

        arc_part = phi_cut-phi_cut_prev
        num_steps = int( math.ceil( l*arc_part/max_step_size ) )
        step_size = arc_part/float( num_steps )

        for j in xrange( num_steps+1 ):
            phi = phi_cut_prev+( float( j )*step_size )
            self.vector_distributions[ a ].add( x0+( l*math.cos( phi ) ), y0+( l*math.sin( phi ) ) )

        return True

    def GetLastMCTruth( self ):
        return self.last_mc_truth

    def Draw( self, path, event_count ):
        if not self.initialized:
            return False

        if type( self.last_mc_truth ) != tuple:
            return False

        graphs = []

        self.update_list_c_i( self.last_mc_truth[ 0 ][ 0 ], self.last_mc_truth[ 0 ][ 1 ], self.last_mc_truth[ 0 ][ 2 ], True )

        for vd in self.vector_distributions:
            if vd.get_size() > 0:
                graphs.append( TGraph( vd.get_size(), vd.get_array_x(), vd.get_array_y() ) )

            else:
                graphs.append( TGraph() )

        limits_x = [ -7.0, 7.0 ]
        limits_y = [ -5.0, 5.0 ]
        filename = path+'/'+self.__class__.__name__+'_arc_'+str( event_count )+'.png'
        draw_graphs( graphs, [ 2+a for a in self.loop_channels ], 0.5, 'x', 'y', limits_x, limits_y, filename )

        for graph in graphs:
            graph.Delete()

        for vd in self.vector_distributions:
            vd.reset()

        hists = []

        for a in self.loop_channels:
            hists.append( TH1F( create_name_hist()              ,
                                ''                              ,
                                2*self.num_bins_margin_t        ,
                                float( -self.num_bins_margin_t ),
                                float( self.num_bins_margin_t  ) ) )

            for n in xrange( 2*self.num_bins_margin_t ):
                hists[ -1 ].SetBinContent( n+1, self.S[ a ][ self.n_trig-self.num_bins_margin_t+n ] )

        filename = path+'/'+self.__class__.__name__+'_trig_'+str( event_count )+'.png'
        draw_hists( hists, [ 2+a for a in self.loop_channels ], 'Time (bin)', '', filename )

        for hist in hists:
            hist.Delete()

        return True
