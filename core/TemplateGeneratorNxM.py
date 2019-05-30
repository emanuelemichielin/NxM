import math
import numpy as np
from scipy.signal import find_peaks, wiener, peak_widths
from matplotlib import pyplot as plt

from ROOT import TGraph, TH1F, TLine

from utils import create_name_hist, get_angle_std, check_angle, vector_distribution, draw_hists, draw_graphs

class TemplateGeneratorNxM:

    def __init__( self, V, calc_r, calc_theta, E_min, E_max, vec_r_lim, mat_theta_lim ):

        # V is the cross-correlation function of the noise
        # calc_r and calc_theta are functors that calculate the radial and azimuthal partitions given the channel amplitudes
        # E_min and E_max define the bin in the total energy to which the templates correspond
        # Therefore, only events with total energy between E_min and E_max will be used to calculate the templates
        # vec_r_lim is a list that contains the limits of the bins in the radial partition
        # mat_theta_lim is a list of lists that, for each bin in the radial partition, contains the limits of the bins in the azimuthal partition
        # Note that, therefore, it is possible to define a segmentation such that the limits of the bins in the azimuthal partition depend on the
        # radial partition
        # As an example, see Fig. 12 of http://titus.stanford.edu/cdms_restricted/elias/2018_12_06-nxm_filter/index.html

        self.num_channels = len( V )
        self.num_bins_t = len( V[ 0 ][ 0 ] )

        self.loop_channels = range( self.num_channels )
        self.loop_bins_t = range( self.num_bins_t )

        self.J = [ [ V[ a ][ a ][ n ] for n in self.loop_bins_t ] for a in self.loop_channels ]

        self.calc_r = calc_r
        self.calc_theta = calc_theta

        self.E_min = E_min
        self.E_max = E_max

        # map_bins_part is the structure that, ultimately, relates each bin in the partition space with the corresponding template
        # It also holds the information about the bin limits provided by vec_r_lim and mat_theta_lim
        # Each element of map_bins_part corresponds to one bin in the radial partition
        # In turn, each such element contains a list where each element corresponds to one bin in the azimuthal partition
        # It is important to notice that the total number of bins registered in map_bins_part is larger than that defined by vec_r_lim and
        # mat_theta_lim
        # The reason is that, in order to handle general segmentations, each bin in the radial partition needs to be complemented by two
        # (or one, in some special cases) further bins that 1) are adjacent to it, and 2) have the same limits in the azimuthal partition as
        # it, called "sideband" bins
        # Note that these sideband bins overlap with other original bins, and therefore one event will, in general, contribute to built
        # templates in more than one bin

        self.map_bins_part = []
        self.cumul_S_fft = []
        self.event_counts = []

        for i in range( len( mat_theta_lim ) ):
            if i > 0:

                # This is the "bottom" sideband bin

                self.map_bins_part.append( ( vec_r_lim[ i-1 ], vec_r_lim[ i ], [] ) )

                for j in range( len( mat_theta_lim[ i ] ) ):
                    self.map_bins_part[ -1 ][ 2 ].append( ( mat_theta_lim[ i ][ j-1 ], mat_theta_lim[ i ][ j ], len( self.cumul_S_fft ) ) )
                    self.cumul_S_fft.append( [ [ 0.0 for n in self.loop_bins_t ] for a in self.loop_channels ] )
                    self.event_counts.append( 0 )

            self.map_bins_part.append( ( vec_r_lim[ i ], vec_r_lim[ i+1 ], [] ) )

            # This is the original bin

            for j in range( len( mat_theta_lim[ i ] ) ):
                self.map_bins_part[ -1 ][ 2 ].append( ( mat_theta_lim[ i ][ j-1 ], mat_theta_lim[ i ][ j ], len( self.cumul_S_fft ) ) )
                self.cumul_S_fft.append( [ [ 0.0 for n in self.loop_bins_t ] for a in self.loop_channels ] )
                self.event_counts.append( 0 )

            if i < len( mat_theta_lim )-1:

                # This is the "top" sideband bin

                self.map_bins_part.append( ( vec_r_lim[ i+1 ], vec_r_lim[ i+2 ], [] ) )

                for j in range( len( mat_theta_lim[ i ] ) ):
                    self.map_bins_part[ -1 ][ 2 ].append( ( mat_theta_lim[ i ][ j-1 ], mat_theta_lim[ i ][ j ], len( self.cumul_S_fft ) ) )
                    self.cumul_S_fft.append( [ [ 0.0 for n in self.loop_bins_t ] for a in self.loop_channels ] )
                    self.event_counts.append( 0 )

        self.vector_distributions = [ vector_distribution() for i in range( len( self.cumul_S_fft )+1 ) ]

    def IncludeEvent( self, S ):
        S_fft = []
        amps = []

        dataS = np.asarray(S)
        dataS = np.sum(dataS, axis=0)
        peaks, properties = find_peaks(dataS.transpose(), prominence = 1, width=10)
        width_half = peak_widths(dataS.transpose(), peaks, rel_height=0.5)
        if (np.mean(dataS[0:3000])>np.mean(dataS[3000:10000]) and np.mean(dataS[0:3000])>np.mean(dataS[10000:15000]) and np.mean(dataS[0:3000])>np.mean(dataS[30000:32000])): return False
        if (np.mean(dataS[0:3000]) > 1.05*np.mean(dataS[30000:32000])) : return False 
        if (1.01*np.mean(dataS[0:3000]) < np.mean(dataS[30000:32000])) : return False

        if ( len(peaks) == 0 or len(peaks) > 1 or any(width_half)>2000) : return False

        #plt.plot(dataS.transpose())
        #plt.plot(peaks, dataS[peaks].transpose(), "xr")
        #plt.show()

        
        tmp_S = np.array(S)
        avg =  np.mean(tmp_S[:,2000:15000],axis=1) 
 
        for i in range(len(S)):
            S[i] = S[i]-avg[i]
        
        for a in self.loop_channels:
            S_fft.append( np.fft.fft( S[ a ] ).tolist() )
            #coefs_fft = [ ( 1.0-( self.J[ a ][ n ]/( np.absolute( S_fft[ -1 ][ n ] )**2 ) ) )*S_fft[ -1 ][ n ] for n in self.loop_bins_t ]
            #amps.append( float( sum( np.fft.ifft( coefs_fft ).tolist() ).real ) )

 
        #using scipy wiener filter
        noise_w = np.array([self.J[a][a] for a in range(len(S))])
        amps = [ sum( wiener( S[ a ], mysize=75, noise = noise_w[a].real ) ) for a in range(len(S)) ]

        E = sum( amps )

        if E < self.E_min or E > self.E_max:
            return False

        r = self.calc_r( amps )
        theta = get_angle_std( self.calc_theta( amps ) )

        found = False

        for i in range( len( self.map_bins_part ) ):
            if r > self.map_bins_part[ i ][ 0 ] and r < self.map_bins_part[ i ][ 1 ]:
                for j in range( len( self.map_bins_part[ i ][ 2 ] ) ):
                    if check_angle( theta, self.map_bins_part[ i ][ 2 ][ j ][ 0 ], self.map_bins_part[ i ][ 2 ][ j ][ 1 ] ):
                        k = self.map_bins_part[ i ][ 2 ][ j ][ 2 ]

                        for a in self.loop_channels:
                            for n in self.loop_bins_t:
                                self.cumul_S_fft[ k ][ a ][ n ] += S_fft[ a ][ n ]

                        self.event_counts[ k ] += 1
                        found = True

                        self.vector_distributions[ k ].add( theta, r )

        self.vector_distributions[ -1 ].add( theta, r )

        return found

    def GetTemplates( self ):
        templates = []

        for i in range( len( self.cumul_S_fft ) ):
            if self.event_counts[ i ] > 0:
                templates.append( [] )

                for a in self.loop_channels:
                    templates[ -1 ].append( [ float( coef.real ) for coef in np.fft.ifft( self.cumul_S_fft[ i ][ a ] ).tolist() ] )

                    for n in self.loop_bins_t:
                        templates[ -1 ][ -1 ][ n ] /= float( self.event_counts[ i ] )

        return templates

    def GetMapBinsPart( self ):
        map_bins_part = []
        template_count = 0

        for i in range( len( self.map_bins_part ) ):
            map_bins_part.append( ( self.map_bins_part[ i ][ 0 ],
                                    self.map_bins_part[ i ][ 1 ],
                                    []                           ) )

            for j in range( len( self.map_bins_part[ i ][ 2 ] ) ):
                if self.event_counts[ self.map_bins_part[ i ][ 2 ][ j ][ 2 ] ] > 0:
                    map_bins_part[ -1 ][ 2 ].append( ( self.map_bins_part[ i ][ 2 ][ j ][ 0 ],
                                                       self.map_bins_part[ i ][ 2 ][ j ][ 1 ],
                                                       template_count                         ) )

                    template_count += 1

                else:

                    # If a bin of the partition space does not has statistics to build the template, then label it with a -1

                    map_bins_part[ -1 ][ 2 ].append( ( self.map_bins_part[ i ][ 2 ][ j ][ 0 ],
                                                       self.map_bins_part[ i ][ 2 ][ j ][ 1 ],
                                                       -1                                     ) )

        return map_bins_part

    def Draw( self, path ):
        graphs = []

        for vd in self.vector_distributions:
            if vd.get_size() > 0:
                graphs.append( TGraph( vd.get_size(), vd.get_array_x(), vd.get_array_y() ) )

            else:
                graphs.append( TGraph() )

        limits_x = [ -math.pi                    , math.pi                       ]
        limits_y = [ self.map_bins_part[ 0 ][ 0 ], self.map_bins_part[ -1 ][ 1 ] ]

        lines = []

        for i in range( len( self.map_bins_part ) ):
            if i%3 == 0:
                r_bot = self.map_bins_part[ i ][ 0 ]
                r_top = self.map_bins_part[ i ][ 1 ]

                for j in range( len( self.map_bins_part[ i ][ 2 ] ) ):
                    theta_clockw = self.map_bins_part[ i ][ 2 ][ j ][ 0 ]
                    theta_anticw = self.map_bins_part[ i ][ 2 ][ j ][ 1 ]

                    lines.append( TLine( theta_clockw, r_bot, theta_clockw, r_top ) )
                    lines.append( TLine( theta_anticw, r_bot, theta_anticw, r_top ) )

                    if theta_anticw > theta_clockw:
                        lines.append( TLine( theta_clockw, r_bot, theta_anticw, r_bot ) )
                        lines.append( TLine( theta_clockw, r_top, theta_anticw, r_top ) )

                    else:
                        lines.append( TLine( theta_clockw, r_bot, math.pi     , r_bot ) )
                        lines.append( TLine( -math.pi    , r_bot, theta_anticw, r_bot ) )
                        lines.append( TLine( theta_clockw, r_top, math.pi     , r_top ) )
                        lines.append( TLine( -math.pi    , r_top, theta_anticw, r_top ) )

        for line in lines:
            line.SetLineColor( 15 )

        for i in range( len( self.map_bins_part ) ):
            if i%3 == 0:
                for j in range( len( self.map_bins_part[ i ][ 2 ] ) ):
                    graphs_subset = [ graphs[ -1 ], graphs[ self.map_bins_part[ i ][ 2 ][ j ][ 2 ] ] ]
                    colors = [ 3, 2 ]

                    if len( self.map_bins_part[ i ][ 2 ] ) > 1:
                        graphs_subset.append( graphs[ self.map_bins_part[ i ][ 2 ][ j-1 ][ 2 ] ] )
                        colors.append( 4 )

                        if len( self.map_bins_part[ i ][ 2 ] ) > 2:
                            graphs_subset.append( graphs[ self.map_bins_part[ i ][ 2 ][ ( j+1 )%len( self.map_bins_part[ i ][ 2 ] ) ][ 2 ] ] )
                            colors.append( 4 )

                    if i > 0:
                        graphs_subset.append( graphs[ self.map_bins_part[ i-1 ][ 2 ][ j ][ 2 ] ] )
                        colors.append( 4 )

                    if i < len( self.map_bins_part )-1:
                        graphs_subset.append( graphs[ self.map_bins_part[ i+1 ][ 2 ][ j ][ 2 ] ] )
                        colors.append( 4 )

                    filename = path+'/'+self.__class__.__name__+'_map_'+str( ( i/3 )+1 )+'_'+str( j+1 )+'.png'
                    draw_graphs( graphs_subset, colors, 0.5, '#theta', 'r', limits_x, limits_y, filename, lines )

        for graph in graphs:
            graph.Delete()

        for vd in self.vector_distributions:
            vd.reset()

        templates = self.GetTemplates()

        hists = []

        for a in self.loop_channels:
            hists.append( TH1F( create_name_hist(), '', self.num_bins_t, 0.0, float( self.num_bins_t ) ) )

        map_bins_part = self.GetMapBinsPart()

        for i in range( len( map_bins_part ) ):
            if i%3 == 0:
                for j in range( len( map_bins_part[ i ][ 2 ] ) ):
                    k = map_bins_part[ i ][ 2 ][ j ][ 2 ]

                    if k > -1:
                        for a in self.loop_channels:
                            for n in self.loop_bins_t:
                                hists[ a ].SetBinContent( n+1, templates[ k ][ a ][ n ] )

                        filename = path+'/'+self.__class__.__name__+'_pulse_'+str( ( i/3 )+1 )+'_'+str( j+1 )+'.png'
                        draw_hists( hists, [ 2+a for a in self.loop_channels ], 'Time (bin)', '', filename, '', legend = True)
                    
                        for hist in hists:
                            hist.Reset()

                    if i > 0:
                        k = map_bins_part[ i-1 ][ 2 ][ j ][ 2 ]

                        if k > -1:
                            for a in self.loop_channels:
                                for n in self.loop_bins_t:
                                    hists[ a ].SetBinContent( n+1, templates[ k ][ a ][ n ] )

                            filename = path+'/'+self.__class__.__name__+'_pulse_'+str( ( i/3 )+1 )+'_'+str( j+1 )+'_bot.png'
                            #draw_hists( hists, [ 2+a for a in self.loop_channels ], 'Time (bin)', '', filename, '')
                    
                            for hist in hists:
                                hist.Reset()

                    if i < len( map_bins_part )-1:
                        k = map_bins_part[ i+1 ][ 2 ][ j ][ 2 ]

                        if k > -1:
                            for a in self.loop_channels:
                                for n in self.loop_bins_t:
                                    hists[ a ].SetBinContent( n+1, templates[ k ][ a ][ n ] )

                            filename = path+'/'+self.__class__.__name__+'_pulse_'+str( ( i/3 )+1 )+'_'+str( j+1 )+'_top.png'
                            #draw_hists( hists, [ 2+a for a in self.loop_channels ], 'Time (bin)', '', filename, '')

                            for hist in hists:
                                hist.Reset()

        for hist in hists:
            hist.Delete()

        return True
