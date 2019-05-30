import math

from ROOT import TLine

def estimate_bins_part( stat, vd ):
    lines = []

    stat_merged = 6*stat
    half_stat_merged = stat_merged/2

    if vd.get_size() > half_stat_merged:
        points = []

        for i in range( vd.get_size() ):
            r = vd.get_array_y()[ i ]
            theta = vd.get_array_x()[ i ]

            if theta < -2.0*math.pi/3.0:
                points.append( [ r, -( theta+( 2.0*math.pi/3.0 ) ) ] )

            elif theta < -math.pi/3.0:
                points.append( [ r, theta+( 2.0*math.pi/3.0 ) ] )

            elif theta < 0.0:
                points.append( [ r, -theta ] )

            elif theta > 2.0*math.pi/3.0:
                points.append( [ r, theta-( 2.0*math.pi/3.0 ) ] )

            elif theta > math.pi/3.0:
                points.append( [ r, -( theta-( 2.0*math.pi/3.0 ) ) ] )

            else:
                points.append( [ r, theta ] )

        vec_r_lim = []
        mat_theta_lim = []

        points.sort( key = lambda point: point[ 0 ] )

        r_min = points[ 0  ][ 0 ]
        r_max = points[ -1 ][ 0 ]

        delta_r = r_max-r_min

        vec_r_lim.append( max( r_min-( 0.1*delta_r ), 0.0 ) )

        prev_point_index = 0
        point_index = half_stat_merged

        i = 0

        while True:
            mat_theta_lim.append( [] )

            if i > 0:
                points_subset = sorted( points[ int(prev_point_index):int(point_index) ], key = lambda point: point[ 1 ] )

                for j in range( 1, i+1 ):
                    if j == 1:
                        mat_theta_lim[ -1 ].append( points_subset[ int(half_stat_merged) ][ 1 ] )

                    else:
                        mat_theta_lim[ -1 ].append( points_subset[ stat_merged*j ][ 1 ] )

            i += 1

            prev_point_index = point_index
            point_index += half_stat_merged+( stat_merged*i )

            if point_index < len( points ):
                vec_r_lim.append( points[ int(prev_point_index) ][ 0 ] )

            else:
                break

        vec_r_lim.append( r_max+( 0.1*delta_r ) )


        for i in range( len( mat_theta_lim ) ):
            sym_theta_lims = []

            for j in range( len( mat_theta_lim[ i ] ) ):
                sym_theta_lims.append( -( mat_theta_lim[ i ][ j ]+( 2.0*math.pi/3.0 ) ) )
                sym_theta_lims.append( mat_theta_lim[ i ][ j ]-( 2.0*math.pi/3.0 )      )
                sym_theta_lims.append( -mat_theta_lim[ i ][ j ]                         )
                sym_theta_lims.append( mat_theta_lim[ i ][ j ]+( 2.0*math.pi/3.0 )      )
                sym_theta_lims.append( -( mat_theta_lim[ i ][ j ]-( 2.0*math.pi/3.0 ) ) )

            for theta in sym_theta_lims:
                mat_theta_lim[ i ].append( theta )

            mat_theta_lim[ i ].append( -math.pi/3.0 )
            mat_theta_lim[ i ].append( math.pi/3.0  )
            mat_theta_lim[ i ].append( math.pi      )

            mat_theta_lim[ i ].sort()
                    
        for i in range( len( mat_theta_lim ) ):
            r_bot = vec_r_lim[ i   ]
            r_top = vec_r_lim[ i+1 ]
 
            for j in range( len( mat_theta_lim[ i ] ) ):
                theta_clockw = mat_theta_lim[ i ][ j-1 ]
                theta_anticw = mat_theta_lim[ i ][ j   ]
 
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



    return lines
