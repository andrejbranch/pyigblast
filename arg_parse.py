import argparse
import textwrap
from Bio import SeqIO
import os
import shutil
from multiprocessing import cpu_count

class blastargument_parser():
  # call on all our file type parsers in the sequence_anlysis_method
    def __init__(self):
        """A customized argument parser that does a LOT of error checking"""
        self.parser = argparse.ArgumentParser(prog="igblast", formatter_class=argparse.RawTextHelpFormatter, description=textwrap.dedent('''\
                                    PyIgBlast
                    __________________________________________________________________________________________\n
                    PyIgBlast calls upon igblastn for nucleotides. Uses multiproecessing to split up fasta file.
                    Parses the output to a csv/tsv/JSON and allows upload to MongoDB or MySQL databases
                    author - Joran Willis
                    '''))

        #Necessary Arguments
        neces = self.parser.add_argument_group(
            title='Necessary', description="These have to be included")
        #query
        neces.add_argument(
            "-q","--query", metavar="query.fasta",required=True,type=self._check_if_fasta,help="The fasta file to be input into igBlast")
        #database path
        neces.add_argument(
            "-d","--db_path",required=True,type=self._check_if_db_exists,help="The database path to the germline repertoire")
        #internal_data path
        neces.add_argument(
            "-i","--internal_data",required=True,type=self._check_if_db_exists,help="The database path to internal data repertoire")

        #recommended options
        recommended = self.parser.add_argument_group(
            title="\nRecommended",description="Not necessary to run but recommended")
        recommended.add_argument(
            "-a","--aux_path",type=self._check_if_aux_path_exists,help="The auxilariay path that contains the frame origins of the germline genes for each repertoire, \
            helps produce translation and other metrics")

        #IGBlast Specif Options
        igspec = self.parser.add_argument_group(
            title="\nIgBlast Sprecific",description="IgBlast Specific Options with a Default")
        igspec.add_argument(
            "-or","--organism",default="human",choices=["human","mouse"],help="The organism repeortire to blast against")
        igspec.add_argument("-nV","--num_v",default=3,type=int,help="How many V-genes to match?")
        igspec.add_argument("-nD","--num_d",default=3,type=int,help="How many D-genes to match?")
        igspec.add_argument("-nJ","--num_j",default=3,type=int,help="How many J-genes to match?")
        igspec.add_argument("-dgm","--d_gene_matches",default=5,type=int,help="How many nuclodtieds in the D-gene must match to call it a hit")
        igspec.add_argument("-s","--domain",default="imgt",choices=["imgt","kabat"],help="Which classification system do you want")
        igspec.add_argument("-sT","--show_translation",default=False,action="store_true",help="Do you want to show the translation of the alignments")

        #General Blast Settings
        general = self.parser.add_argument_group(
            title="\nGeneral Settings",description="General Settings for Blast")
        general.add_argument("-x",'--executable',type=self._check_if_executable_exists,help="The location of the executable, default is /usr/bin/igblastn")
        general.add_argument("-o","--out",help="output file prefix",default="igblast_out_")
        general.add_argument("-e","--e_value",type=str,default="1e-15",help="Real value for excpectation value threshold in blast, put in scientific notation")
        general.add_argument("-w","--word_size",type=int,default=4,help="Word size for wordfinder algorithm")
        general.add_argument("-pm","--penalty_mismatch",type=int,default=0,help="Penalty for nucleotide mismatch")
        general.add_argument("-rm","--reward_match",type=int,default=0,help="Reward for nucleotide match")
        general.add_argument("-mT","--max_target_seqs",type=int,default=500,help="Maximum number of alingned sequences to iterate through at a time")
        general.add_argument("-nP","--num_procs",type=int,default=cpu_count(),help="How many do you want to split the job across, default is the number of processors")

        formatter = self.parser.add_argument_group(
            title="Formatting Options",description="Formatting options mostly available"
            )
        formatter.add_argument("-f","--format_options",type=str,default="default",help="default is a tab seperated format of\n\n\
         qseqid sseqid pident length mismatch gapopen qstart qend sstart send\n\n\
         The format file is in the database path as format_template.txt. Uncomment out the metrics you want to use")


        #one special boolean case
        self.show_translation = False
        #return the arguments
        self.args = self.parser.parse_args()
        #get them ready to ship out
        self._make_args_dict()


    #helper functions
    def _check_if_fasta(self,f_file):
        try:
            SeqIO.parse(f_file, "fasta").next()
            return f_file

            #return SeqIO.parse(f_file,"fasta")
        except StopIteration:
            msg = "{0} is not a fasta file\n".format(f_file)
            raise argparse.ArgumentTypeError(msg)

    def _check_if_executable_exists(self,x_path):
        if not os.path.exists(x_path):
            msg = "path to executable {0} does not exist, use -h for help\n".format(x_path)
            raise argparse.ArgumentTypeError(msg)
        if not os.access(x_path, os.R_OK):
            msg1 = "executable {0} does not permission to run\n".format(x_path)
            raise argparse.ArgumentTypeError(msg1)
        else:
            return x_path


    def _check_if_db_exists(self,db_path):
        if os.path.exists(db_path):
            return db_path
        else:
            msg = "{0} path for does not exist for database\n".format(db_path)
            raise argparse.ArgumentTypeError(msg)

    def _check_if_aux_path_exists(self,aux_path):
        if os.path.exists(aux_path):
            return aux_path
        else:
            msg = "{0} path for aux files does not exist\n".format(aux_path)
            raise argparse.ArgumentTypeError(msg)

    def _make_args_dict(self):
        #copy internal data directory to current location
        #shutil.copytree(self.args.internal_data,'.')
        try:
            shutil.copytree(self.args.internal_data,'./internal_data')
        except OSError:
            print "Internal Data direcotry file exists in this directory, skipping..."


        self.args_dict = {
                      '-query':self.args.query,
                      '-organism':self.args.organism,
                      '-num_alignments_V':self.args.num_v,
                      '-num_alignments_D':self.args.num_d,
                      '-num_alignments_J':self.args.num_j,
                      '-min_D_match':self.args.d_gene_matches,
                      '-domain_system':self.args.domain,
                      '-out':self.args.out,
                      '-evalue':self.args.e_value,
                      '-word_size':self.args.word_size,
                      '-max_target_seqs':self.args.max_target_seqs,
                      '-germline_db_V':"{0}{1}_gl_V".format(self.args.db_path,self.args.organism),
                      '-germline_db_D':"{0}{1}_gl_D".format(self.args.db_path,self.args.organism),
                      '-germline_db_J':"{0}{1}_gl_J".format(self.args.db_path,self.args.organism)

        }
        #add bool opition
        if self.args.penalty_mismatch:
            self.args_dict['-penalty'] = self.args.penalty_mismatch
        if self.args.reward_match:
            self.args_dict['-reward'] = self.args.reward_match
        if self.args.show_translation:
            self.show_translation = True
        if self.args.aux_path:
            self.args_dict['-auxiliary_data'] = "{0}{1}_gl.aux".format(self.args.aux_path,self.args.organism)

        #add formatting option
        if self.args.format_options == 'default':
            self.args_dict['-outfmt'] = 7
        else:
            self.args.format_options = self._check_if_db_exists(self.args.format_options)
            formatting_titles = []
            for line in open(self.args.format_options).readlines():
                if line.startswith("#"):
                    continue
                else:
                    formatting_titles.append(line.split()[0])

            format = "7 " + " ".join(formatting_titles)
            self.args_dict['-outfmt'] = format


    #only non memeber functions needed
    def return_parsed_args(self):
        return self.args_dict

    def return_command_line(self):
        '''return solid strin of command line'''
        return ' '.join(self.return_command_line_from_dict(self.args_dict))

    def return_command_line_from_dict(self,cline_dict):
        '''return command line as a list to put in subprocess
        --args
        cline_dict - The command line dictionary to return. We add in the executable'''
        if self.args.executable:
            cline = [self.args.executable]
        else:
            cline = [self._check_if_executable_exists("/usr/bin/igblastn")]
        for command in cline_dict:
            cline.append(str(command))
            cline.append(str(self.args_dict[command]))
        return cline


if __name__ == '__main__':
    args = blastargument_parser().return_command_line()
    print args
