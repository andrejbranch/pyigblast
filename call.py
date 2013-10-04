import subprocess
import sys
import os
import glob

class call_igblast():
    '''A class for handling the input to igblasn and calling on the executable

    ---args 
    input - The fastafile of the input 
    d_locatin - The location of the database folder 
    --KeyWord
    location - The location of the igBlastExecutable (default location="/usr/bin/igblastn")'''
    
    def __init__(self,input,d_location,location="/usr/bin/igblastn"):
        #all I need now is the location of the executable
        self.location = location
        self.database_folder = d_location
        if os.path.exists(self.database_folder):
            self._parse_db()
        else:
            print 'Check location of database folder {0}'.format(input)
        if os.path.exists(input):
            self.input = input
        else:
            print 'Check location of input file {0}'.format(input)
            sys.exit()
    def call_process(self):
        pass

    def _parse_db(self,):
        files = glob.glob(self.database_folder)
        for file in files:
            print file
