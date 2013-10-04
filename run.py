import sys

from call import call_igblast as cig 
input_file  = sys.argv[1]
database_location = sys.argv[2]
aux_location = sys.argv[3]

run_ig = cig(input_file,database_location,aux_location)
run_ig.run()