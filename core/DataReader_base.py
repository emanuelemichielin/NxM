class DataReader_base:

    def __init__( self ):
        self.S = None

    def OpenFile( self, filename ):

        # This function needs to be overloaded by the user in order to open the data files of interest

        return False

    def LoadEvent( self ):

        # This function needs to be overloaded by the user in order to read events from the data files of interest

        return False

    def GetTraces( self ):
        return self.S

    def CloseFile( self ):

        # This function needs to be overloaded by the user in order toclose the data files of interest

        self.S = None

        return True
