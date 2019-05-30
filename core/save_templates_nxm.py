import array

def save_templates_nxm( templates, E_min, E_max, map_bins_part, filename ):
    filepointer = open( filename, 'wb')

    num_templates = len( templates )
    array_num_templates = array.array( 'l', [ num_templates ] )
    array_num_templates.tofile( filepointer )

    num_channels = len( templates[ 0 ] )
    array_num_channels = array.array( 'l', [ num_channels ] )
    array_num_channels.tofile( filepointer )

    num_bins_t = len( templates[ 0 ][ 0 ] )
    array_num_bins_t = array.array( 'l', [ num_bins_t ] )
    array_num_bins_t.tofile( filepointer )

    for i in range( num_templates ):
        for a in range( num_channels ):
            array_template = array.array( 'd', templates[ i ][ a ] )
            array_template.tofile( filepointer )

    array_E_min = array.array( 'd', [ E_min ] )
    array_E_min.tofile( filepointer )

    array_E_max = array.array( 'd', [ E_max ] )
    array_E_max.tofile( filepointer )

    num_bins_r = len( map_bins_part )
    array_num_bins_r = array.array( 'l', [ num_bins_r ] )
    array_num_bins_r.tofile( filepointer )

    for i in range( num_bins_r ):
        array_r_bot = array.array( 'd', [ map_bins_part[ i ][ 0 ] ] )
        array_r_bot.tofile( filepointer )

        array_r_top = array.array( 'd', [ map_bins_part[ i ][ 1 ] ] )
        array_r_top.tofile( filepointer )

        num_bins_theta = len( map_bins_part[ i ][ 2 ] )
        array_num_bins_theta = array.array( 'l', [ num_bins_theta ] )
        array_num_bins_theta.tofile( filepointer )

        for j in range( num_bins_theta ):
            array_theta_clockw = array.array( 'd', [ map_bins_part[ i ][ 2 ][ j ][ 0 ] ] )
            array_theta_clockw.tofile( filepointer )

            array_theta_anticw = array.array( 'd', [ map_bins_part[ i ][ 2 ][ j ][ 1 ] ] )
            array_theta_anticw.tofile( filepointer )

            array_template_index = array.array( 'l', [ map_bins_part[ i ][ 2 ][ j ][ 2 ] ] )
            array_template_index.tofile( filepointer )

    filepointer.close()
    return True
