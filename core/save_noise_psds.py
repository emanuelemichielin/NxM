import array
import numpy as np

def save_noise_psds( noise_psds, filename ):
    filepointer = open( filename, 'wb')

    num_channels = len( noise_psds )
    array_num_channels = array.array( 'l', [ num_channels ] )
    array_num_channels.tofile( filepointer )

    num_bins_t = len( noise_psds[ 0 ][ 0 ] )
    array_num_bins_t = array.array( 'l', [ num_bins_t ] )
    array_num_bins_t.tofile( filepointer )

    for a in range( num_channels ):
        for b in range( num_channels ):
            array_psds_re = array.array( 'd', [ noise_psds[ a ][ b ][ n ].real for n in range( num_bins_t ) ] )
            array_psds_re.tofile( filepointer )

            array_psds_im = array.array( 'd', [ noise_psds[ a ][ b ][ n ].imag for n in range( num_bins_t ) ] )
            array_psds_im.tofile( filepointer )

    filepointer.close()
    return True
