import array
import math

from ROOT import gDirectory, gPad, gROOT, gStyle, TCanvas, THStack, TTree, TLegend

def create_name_hist():
    return 'h'+str( gDirectory.GetList().GetEntries()+1 )

def get_angle_std( angle ):

    # For a given input angle, calculates the equivalent angle in the range (-pi, pi]

    angle_std = angle

    while angle_std > math.pi:
        angle_std -= 2.0*math.pi

    while angle_std <= -math.pi:
        angle_std += 2.0*math.pi

    return angle_std

def check_angle( theta, limit_theta_clockw, limit_theta_anticw ):

    # Checks whether a given angle is contained within a given interval
    # The interval is defined by providing the limits that are situated clockwise/anticlockwise from its center

    if limit_theta_anticw > limit_theta_clockw:
        if theta > limit_theta_clockw and theta <= limit_theta_anticw:
            return True

    else:
        if theta > limit_theta_clockw or theta <= limit_theta_anticw:
            return True

    return False

def calculate_delta_angle( theta_clockw, theta_anticw ):
    if theta_anticw > theta_clockw:
        return theta_anticw-theta_clockw

    else:
        return ( 2.0*math.pi )-( theta_clockw-theta_anticw )

def interpolate_parab( y1, y2, y3 ):
    return [ ( 0.5*( y1+y3 ) )-y2,
             0.5*( y3-y1 )       ,
             y2                   ]

def interpolate_parab_xy( x, y, constrained = False ):
    if constrained:
        if len( x ) != 2 or len( y ) != 2:
            return None

        coefs = np.solve( [ [ x[ 0 ]*x[ 0 ], 1.0 ],
                           [ x[ 1 ]*x[ 1 ], 1.0 ] ], vec ).tolist() 

        return [ coefs[ 0 ], 0.0, coefs[ 1 ] ]

    else:
        if len( x ) != 3 or len( y ) != 3:
            return None

        return np.solve( [ [ x[ 0 ]*x[ 0 ], x[ 0 ], 1.0 ],
                           [ x[ 1 ]*x[ 1 ], x[ 1 ], 1.0 ],
                           [ x[ 2 ]*x[ 2 ], x[ 2 ], 1.0 ] ], vec ).tolist()

def configure_draw():
    gROOT.SetStyle( 'Plain' )

    gStyle.SetEndErrorSize( 0. )
    gStyle.SetOptStat( 0 )
    gStyle.SetPalette( 1 )
    gStyle.SetTitleBorderSize( 0 )
    gStyle.SetTitleFont( 22, '' )
    gStyle.SetTitleSize( 0.06, '' )
    gStyle.SetTitleX( 0.1 )
    gStyle.SetTitleW( 0.8 )

class vector_distribution:

    def __init__( self ):
        self.list_x = []
        self.list_y = []

    def add( self, x, y ):
        self.list_x.append( x )
        self.list_y.append( y )

    def get_size( self ):
        return len( self.list_x )

    def get_array_x( self ):
        return array.array( 'd', self.list_x )

    def get_array_y( self ):
        return array.array( 'd', self.list_y )

    def reset( self ):
        self.list_x = []
        self.list_y = []

def draw_hists( hists, colors, title_x, title_y, filename, a = '', log_y = False, legend = False ):
    if len( hists ) == 0 or len( colors ) != len( hists ):
        return False

    for ind in range( len( hists ) ):
        hists[ ind ].SetLineColor( colors[ ind ] )
        hists[ ind ].SetMarkerColor( colors[ ind ] )
        
    c = TCanvas( 'c', '',800, 600 )

    if log_y:
        gPad.SetLogy( 1 )

    if len( hists ) == 1:
        hists[ 0 ].GetXaxis().SetLabelFont( 22 )
        hists[ 0 ].GetXaxis().SetLabelOffset( 0.012 )
        hists[ 0 ].GetXaxis().SetLabelSize( 0.035 )
        hists[ 0 ].GetXaxis().SetTitle( title_x )
        hists[ 0 ].GetXaxis().SetTitleFont( 22 )
        hists[ 0 ].GetXaxis().SetTitleOffset( 1.34 )
        hists[ 0 ].GetXaxis().SetTitleSize( 0.04 )

        hists[ 0 ].GetYaxis().SetLabelFont( 22 )
        hists[ 0 ].GetYaxis().SetLabelOffset( 0.012 )
        hists[ 0 ].GetYaxis().SetLabelSize( 0.035 )
        hists[ 0 ].GetYaxis().SetTitle( title_y )
        hists[ 0 ].GetYaxis().SetTitleFont( 22 )
        hists[ 0 ].GetYaxis().SetTitleOffset( 1.4 )
        hists[ 0 ].GetYaxis().SetTitleSize( 0.04 )

        hists[ 0 ].Draw()
        c.SaveAs( filename )

    else:
        leg = TLegend(.55,.5,.9,.88)
        hstack = THStack( 'hstack', '' )

        for hist in hists:
            hstack.Add( hist )
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        for idx, hist in enumerate(hists):
            channel = "Channel "+str(a) + " " + str(idx)
            leg.AddEntry(hist,channel,"l")


        hstack.Draw( 'nostack' )
        c.Update()

        hstack.GetXaxis().SetLabelFont( 22 )
        hstack.GetXaxis().SetLabelOffset( 0.012 )
        hstack.GetXaxis().SetLabelSize( 0.035 )
        hstack.GetXaxis().SetTitle( title_x )
        hstack.GetXaxis().SetTitleFont( 22 )
        hstack.GetXaxis().SetTitleOffset( 1.34 )
        hstack.GetXaxis().SetTitleSize( 0.04 )

        hstack.GetYaxis().SetLabelFont( 22 )
        hstack.GetYaxis().SetLabelOffset( 0.012 )
        hstack.GetYaxis().SetLabelSize( 0.035 )
        hstack.GetYaxis().SetTitle( title_y )
        hstack.GetYaxis().SetTitleFont( 22 )
        hstack.GetYaxis().SetTitleOffset( 1.4 )
        hstack.GetYaxis().SetTitleSize( 0.04 )

        #hstack.SetMinimum( gPad.GetUymin() )
        #hstack.SetMaximum( gPad.GetUymax() )
        hstack.Draw( 'nostack' )
        if (legend): leg.Draw()
        c.SaveAs( str(filename) )
        hstack.Delete()
        leg.Delete()

    if log_y:
        gPad.SetLogy( 0 )
    c.Close()
    return True

def draw_hist_2D( hist, title_x, title_y, filename ):
    c = TCanvas( 'c', '' )

    hist.GetXaxis().SetLabelFont( 22 )
    hist.GetXaxis().SetLabelOffset( 0.012 )
    hist.GetXaxis().SetLabelSize( 0.035 )
    hist.GetXaxis().SetTitle( title_x )
    hist.GetXaxis().SetTitleFont( 22 )
    hist.GetXaxis().SetTitleOffset( 1.34 )
    hist.GetXaxis().SetTitleSize( 0.04 )

    hist.GetYaxis().SetLabelFont( 22 )
    hist.GetYaxis().SetLabelOffset( 0.012 )
    hist.GetYaxis().SetLabelSize( 0.035 )
    hist.GetYaxis().SetTitle( title_y )
    hist.GetYaxis().SetTitleFont( 22 )
    hist.GetYaxis().SetTitleOffset( 1.4 )
    hist.GetYaxis().SetTitleSize( 0.04 )

    hist.Draw( 'colz' )
    c.SaveAs( filename )

    c.Close()
    return True

def draw_graphs( graphs, colors, dotsize, title_x, title_y, limits_x, limits_y, filename, lines = [] ):
    if len( graphs ) == 0 or len( colors ) != len( graphs ) or len( limits_x ) != 2 or len( limits_y ) != 2:
        return False

    for ind in range( len( graphs ) ):
        graphs[ ind ].SetMarkerColor( colors[ ind ] )
        graphs[ ind ].SetMarkerSize( dotsize )
        graphs[ ind ].SetMarkerStyle( 8 )

    c = TCanvas( 'c', '' )

    ind = 0

    while True:
        if graphs[ ind ].GetN() > 0:
            graphs[ ind ].GetXaxis().SetLabelFont( 22 )
            graphs[ ind ].GetXaxis().SetLabelOffset( 0.012 )
            graphs[ ind ].GetXaxis().SetLabelSize( 0.035 )
            graphs[ ind ].GetXaxis().SetTitle( title_x )
            graphs[ ind ].GetXaxis().SetTitleFont( 22 )
            graphs[ ind ].GetXaxis().SetTitleOffset( 1.34 )
            graphs[ ind ].GetXaxis().SetTitleSize( 0.04 )

            graphs[ ind ].GetYaxis().SetLabelFont( 22 )
            graphs[ ind ].GetYaxis().SetLabelOffset( 0.012 )
            graphs[ ind ].GetYaxis().SetLabelSize( 0.035 )
            graphs[ ind ].GetYaxis().SetTitle( title_y )
            graphs[ ind ].GetYaxis().SetTitleFont( 22 )
            graphs[ ind ].GetYaxis().SetTitleOffset( 1.4 )
            graphs[ ind ].GetYaxis().SetTitleSize( 0.04 )

            graphs[ ind ].SetTitle( '' )

            graphs[ ind ].GetXaxis().SetLimits( limits_x[ 0 ], limits_x[ 1 ] )

            graphs[ ind ].SetMinimum( limits_y[ 0 ] )
            graphs[ ind ].SetMaximum( limits_y[ 1 ] )

            graphs[ ind ].Draw( 'ap' )
            break

        ind += 1

    if ind == len( graphs ):
        return False

    for graph in graphs[ :ind ]:
        if graph.GetN() > 0:
            graph.Draw( 'p' )

    for graph in graphs[ ind+1: ]:
        if graph.GetN() > 0:
            graph.Draw( 'p' )

    for line in lines:
        line.Draw()

    c.SaveAs( filename )

    c.Close()
    return True

class tree_manager:

    def __init__( self, name ):
        self.tree = TTree( name, '' )
        self.arrays = {}

    def Branch( self, key ):
        self.arrays[ key ] = array.array( 'd', [ 0.0 ] )
        return self.tree.Branch( key, self.arrays[ key ], key+'/D' )

    def __setitem__( self, key, value ):
        self.arrays[ key ][ 0 ] = value

    def Fill( self ):
        return self.tree.Fill()
    
    def Write( self ):
        return self.tree.Write()
