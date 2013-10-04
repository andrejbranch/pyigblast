#import sys
from call import call_igblast as cig 
from arg_parse import blastargument_parser

ap = blastargument_parser().return_parsed_args()
run_ig = cig(ap.query,ap.db_path)
run_ig.run()