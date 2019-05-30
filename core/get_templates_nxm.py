import array

def get_templates_nxm( filename ):
    templates = []
    map_bins_part = []

    filepointer = open( filename, 'rb' )

    data = [ filepointer.read( 8 ), filepointer.read( 8 ), filepointer.read( 8 ) ]

    if len( data[ 0 ] ) < 8 or len( data[ 1 ] ) < 8 or len( data[ 2 ] ) < 8:
        return None

    array_num_templates = array.array( 'l' )
    array_num_templates.frombytes( data[ 0 ] )
    num_templates = array_num_templates[ 0 ]

    array_num_channels = array.array( 'l' )
    array_num_channels.frombytes( data[ 1 ] )
    num_channels = array_num_channels[ 0 ]

    array_num_bins_t = array.array( 'l' )
    array_num_bins_t.frombytes( data[ 2 ] )
    num_bins_t = array_num_bins_t[ 0 ]

    for i in range( num_templates ):
        templates.append( [] )

        for a in range( num_channels ):
            data = [ filepointer.read( 8*num_bins_t ) ]

            if len( data[ 0 ] ) < 8*num_bins_t:
                return None

            array_template = array.array( 'd' )
            array_template.frombytes( data[ 0 ] )
            templates[ -1 ].append( array_template.tolist() )

    data = [ filepointer.read( 8 ), filepointer.read( 8 ), filepointer.read( 8 ) ]

    if len( data[ 0 ] ) < 8 or len( data[ 1 ] ) < 8 or len( data[ 2 ] ) < 8:
        return None

    array_E_min = array.array( 'd' )
    array_E_min.frombytes( data[ 0 ] )
    E_min = array_E_min[ 0 ]

    array_E_max = array.array( 'd' )
    array_E_max.frombytes( data[ 1 ] )
    E_max = array_E_max[ 0 ]

    array_num_bins_r = array.array( 'l' )
    array_num_bins_r.frombytes( data[ 2 ] )
    num_bins_r = array_num_bins_r[ 0 ]

    for i in range( num_bins_r ):
        data = [ filepointer.read( 8 ), filepointer.read( 8 ), filepointer.read( 8 ) ]

        if len( data[ 0 ] ) < 8 or len( data[ 1 ] ) < 8 or len( data[ 2 ] ) < 8:
            return None

        array_r_bot = array.array( 'd' )
        array_r_bot.frombytes( data[ 0 ] )

        array_r_top = array.array( 'd' )
        array_r_top.frombytes( data[ 1 ] )

        map_bins_part.append( ( array_r_bot[ 0 ], array_r_top[ 0 ], [] ) )

        array_num_bins_theta = array.array( 'l' )
        array_num_bins_theta.frombytes( data[ 2 ] )
        num_bins_theta = array_num_bins_theta[ 0 ]

        for j in range( num_bins_theta ):
            data = [ filepointer.read( 8 ), filepointer.read( 8 ), filepointer.read( 8 ) ]

            if len( data[ 0 ] ) < 8 or len( data[ 1 ] ) < 8 or len( data[ 2 ] ) < 8:
                return None

            array_theta_clockw = array.array( 'd' )
            array_theta_clockw.frombytes( data[ 0 ] )

            array_theta_anticw = array.array( 'd' )
            array_theta_anticw.frombytes( data[ 1 ] )

            array_template_index = array.array( 'l' )
            array_template_index.frombytes( data[ 2 ] )

            map_bins_part[ -1 ][ 2 ].append( ( array_theta_clockw[ 0 ]  ,
                                               array_theta_anticw[ 0 ]  ,
                                               array_template_index[ 0 ] ) )

    if len( filepointer.read() ) > 0:
        return None

    filepointer.close()
    return templates, E_min, E_max, map_bins_part
