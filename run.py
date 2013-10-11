from arg_parse import blastargument_parser
import subprocess as sp
import multiprocessing as mp
import glob
import os
import output_parser
import sys
import gzip

try:
    import Bio
except ImportError("Trouble Installing BioPython:"):
    print "Can't find BioPython Module in this path. PyIgBlast is dependent on Biopython"
from time import time

# setup global class to pass around the command line
arg_parser = blastargument_parser()
# holds the arguments in a dictionary
arg_dict = arg_parser.return_parsed_args()


def split_file(num_procs):
        '''Split the file name by the number of processors you specify

        arguments (num_procs - The amount of processors you are running, lets the function know how much to split up the file into with one file per processor)'''

        file_prefix = arg_dict['-query']
        print "Counting entries in fasta files..."
        parent_file = list(Bio.SeqIO.parse(arg_dict['-query'], 'fasta'))
        length_parent_file = len(parent_file)
        print "{0} in fasta file".format(length_parent_file)

        # can be aproximation
        files_per_tmp = float(length_parent_file) / float(num_procs)
        print "{0} processors, blasting {1} entries per processor".format(num_procs, files_per_tmp)

        # carries all of our records to be written
        joiner = []
        file_counter = 1
        num = 1

        # enumerate huge fasta
        for record in parent_file:
            # append records to our list holder
            joiner.append(">" + record.id + "\n" + str(record.seq))
            # if we have reached the maximum numbers to be in that file, write
            # to a file, and then clear
            if num > files_per_tmp:
                joiner.append("")
                with open(file_prefix + str(file_counter) + ".tmp.fasta", 'w') as f:
                    f.write("\n".join(joiner))
                # change file name,clear record holder, and change the file
                # count
                joiner = []
                file_counter += 1
                num = 1
            else:
                num += 1
        if joiner:
            # for left over fasta entries, very important or else they will
            # just hang out in limbo
            joiner.append("")
            with open(file_prefix + str(file_counter) + ".tmp.fasta", 'w') as f:
                f.write("\n".join(joiner))
        return length_parent_file


def run(file_num):
    '''run, give it a fasta file'''
    # append the file name before running the file
    arg_dict['-query'] = file_num
    arg_dict['-out'] = file_num + ".blast_out"
    cline = arg_parser.return_command_line_from_dict(arg_dict)
    # print " ".join(cline)

    # special case - can't figure how else to add this option
    if arg_parser.show_translation:
        cline.append("-show_translation")
    # call on igBlast
    sp.call(cline)

    # bounty of output options
    zip_bool = arg_parser.args.zip
    json_bool = arg_parser.args.json
    concat_bool = arg_parser.args.concatenate
    blast_out = file_num + ".blast_out"
    zip_out = file_num + "blast_out.gz"
    json_out = arg_parser.args.json_prefix + file_num + ".json"
    if json_bool and zip_bool:
        output_parser.igblast_output(blast_out, json_out, gz=zip_bool)
        if concat_bool:
            os.remove(blast_out)
            os.remove(file_num)
    elif json_bool and concat_bool:
        # wants the json back and to remove the temp
        output_parser.igblast_output(blast_out, json_out, gz=zip_bool)
        os.remove(blast_out)
        os.remove(file_num)
    elif json_bool and not concat_bool:
        # wan't the json but don't want the files removed
        output_parser.igblast_output(blast_out, json_out, gz=zip_bool)
    elif zip_bool and not json_bool:
        #"want the zipped blast outs"
        f_in = open(blast_out, 'rb')
        f_out = gzip.open(zip_out, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        if concat_bool:
            os.remove(file_num)
            os.remove(blast_out)
    elif concat_bool and not json_bool:
        #"don't want to zip up but want to remove"
        os.remove(file_num)


def concat_files():
    '''goes through output formats and concats all files'''
    # concat all files and remove them after
    zip_bool = arg_parser.args.zip
    json_bool = arg_parser.args.json
    concat_bool = arg_parser.args.concatenate
    json_out = arg_parser.args.json_prefix
    global_prefix = arg_dict['-out']
    common_name = arg_dict['-query']
    if zip_bool and json_bool and concat_bool:
        zipped_and_json = glob.glob("*" + common_name + "*" + "json.gz")
        with gzip.open(json_out + ".json.gz", 'wb') as gf:
            for file in zipped_and_json:
                f_in = gzip.open(file, 'rb')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)

    elif json_bool and concat_bool and not zip_bool:
        json_docs = glob.glob("*" + common_name + "*" + "json")
        with open(global_prefix + ".json", 'w') as gf:
            for file in json_docs:
                f_in = open(file, 'r')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)

    elif zip_bool and concat_bool and not json_bool:
        zipped = glob.glob("*" + common_name + "*" + "blast_out.gz")
        with gzip.open(global_prefix + ".blast_out.gz", 'wb') as gf:
            for file in zipped:
                f_in = gzip.open(file, 'rb')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)

    elif concat_bool and not zip_bool and not json_bool:
        blast_out = glob.glob("*" + common_name + "*" + "blast_out")
        with open(global_prefix + ".blast_out", 'w') as gf:
            for file in blast_out:
                f_in = open(file, 'r')
                gf.writelines(f_in)
                f_in.close()
                os.remove(file)


def main():
    now = time()
    num_procs = arg_parser.args.num_procs - 1
    mp_pool = mp.Pool(processes=num_procs)
    entry_num = split_file(num_procs)
    tmp_file = glob.glob(arg_dict['-query'] + "*.tmp.fasta")
    mp_pool.map(run, tmp_file)
    concat_files()
    then = time()
    print "Process of {0} seqeunces using {1} processors took {2} seconds".format(entry_num, num_procs, then - now)
