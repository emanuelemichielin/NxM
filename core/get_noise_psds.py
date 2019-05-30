import array

def get_noise_psds( filename ):
    noise_psds = []

    filepointer = open( filename, 'rb' )

    data = [ filepointer.read( 8 ), filepointer.read( 8 ) ]

    if len( data[ 0 ] ) < 8 or len( data[ 1 ] ) < 8:
        return None

    array_num_channels = array.array( 'l' )
    array_num_channels.frombytes( data[ 0 ] )
    num_channels = array_num_channels[ 0 ]

    array_num_bins_t = array.array( 'l' )
    array_num_bins_t.frombytes( data[ 1 ] )
    num_bins_t = array_num_bins_t[ 0 ]

    for a in range( num_channels ):
        noise_psds.append( [] )

        for b in range( num_channels ):
            data = [ filepointer.read( 8*num_bins_t ), filepointer.read( 8*num_bins_t ) ]

            if len( data[ 0 ] ) < 8*num_bins_t or len( data[ 1 ] ) < 8*num_bins_t:
                return None

            array_psd_re = array.array( 'd' )
            array_psd_re.frombytes( data[ 0 ] )

            array_psd_im = array.array( 'd' )
            array_psd_im.frombytes( data[ 1 ] )

            noise_psds[ -1 ].append( [ array_psd_re[ n ]+( 1.0j*array_psd_im[ n ] ) for n in range( num_bins_t ) ] )

    if len( filepointer.read() ) > 0:
        return None

    filepointer.close()
    return noise_psds
