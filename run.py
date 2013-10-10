from arg_parse import blastargument_parser
import subprocess as sp
import multiprocessing as mp
import glob
import os
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


def remove_temp_files(files_to_remove):
    for file in files_to_remove:
        os.remove(file)


def run(file_num):
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


def concat_files(out_files):
    # concat all files and remove them after
    out_files = [i + ".blast_out" for i in out_files]
    with open(arg_dict['-out'] + "all.blast_out", 'w') as outfile:
        for fname in out_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            os.remove(fname)


def main():
    now = time()
    num_procs = arg_parser.args.num_procs - 1
    mp_pool = mp.Pool(processes=num_procs)
    entry_num = split_file(num_procs)
    tmp_file = glob.glob(arg_dict['-query'] + "*.tmp.fasta")
    mp_pool.map(run, tmp_file)
    remove_temp_files(tmp_file)
    concat_files(tmp_file)
    then = time()
    print "Process of {0} seqeunces using {1} processors took {2} seconds".format(entry_num, num_procs, then - now)
