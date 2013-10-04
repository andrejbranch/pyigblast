import subprocess
import sys
import os
import glob
#import re



class call_igblast():
    '''A class for handling the input to igblasn and calling on the executable

    ---args 
    input - The fastafile of the input 
    d_locatin - The location of the database folder 
    --KeyWord
    location - The location of the igBlastExecutable (default location="/usr/bin/igblastn")'''
    
    def __init__(self,input,d_location,aux_location,location="/usr/bin/igblastn",species='human'):
        #all I need now is the location of the executable
        self.location = location
        self.database_folder = d_location
        self.species = species
        #database input
        if os.path.exists(self.database_folder):
            self._parse_db()
            self._get_database_prepends()
        else:
            print 'Check location of database folder {0}'.format(input)
            sys.exit()

        #query input
        if os.path.exists(input):
            self.input = input
        else:
            print 'Check location of input file {0}'.format(input)
            sys.exit()

        #aux location
        if os.path.exists(aux_location):
            self.aux_location = aux_location
            self._parse_aux()
        else:
            print "Check location of location of auxiliry folder {0}\n".format(aux_location)
            checker = raw_input("Press 1 to continue without aux file\n")
            if checker != "1":
                print "Cancelling Job"
                sys.exit()

    def run(self):

        #set up args
        args = [self.location,'-query',self.input,'-germline_db_V',self.prepends[self.species+"_v"],
        '-germline_db_D',self.prepends[self.species+"_d"],'-germline_db_J',self.prepends[self.species+"_j"],
        "-outfmt","7","-show_translation"
        ]

        #did they specify and aux file?
        if self.aux_file:
            args += ['-auxiliary_data',self.aux_file]
        print "Run:\n{0}".format(" ".join(args))
        
        #call
        subprocess.call(args)


    def _parse_db(self):
        '''parse full database files'''
        self.files = glob.glob(self.database_folder+"/*")
        self.human_database_files = []
        self.mouse_database_files = []
        self.custom_human_database_files = []
        for file in self.files:
            start_name = os.path.basename(file).split("_")[0]
            if start_name == "human":
                self.human_database_files.append(file)
            elif start_name == "mouse":
                self.mouse_database_files.append(file)
            elif start_name == "custom":
                self.custom_human_database_files.append(file)

    def _get_database_prepends(self):
        human_common_prefix = os.path.commonprefix(self.human_database_files)
        mouse_common_prefix = os.path.commonprefix(self.mouse_database_files)
        custom_common_prefix = os.path.commonprefix(self.custom_human_database_files)
        self.prepends = {'human_v':human_common_prefix+"V",
                         'human_d':human_common_prefix+"D",
                         'human_j':human_common_prefix+"J",
                         'mouse_v':mouse_common_prefix+"V",
                         'mouse_d':mouse_common_prefix+"D",
                         'mouse_j':mouse_common_prefix+"J",
                         'custom_v':custom_common_prefix+"V"
        }

    def _parse_aux(self):
        self.aux_file = '{0}{1}_g1_aux'.format(self.aux_location, self.species)