
import os, sys, re, atexit, subprocess, time
from io import StringIO
from Bio.Data import CodonTable
from Bio import SeqIO
from checkSplice import checkSplice
from checkFunctionality import checkFunctionality
from assignNames import assignNames
from parseRSS import parseRSS
from blast2bed import blast2bed, blastOnly
from parseTRA import parseTRA

########## COMMAND LOGGING ############
global printLog
printLog = False #this tells us whether or not to print log info on exit (skip if program was called with -h)


#overload default error handling so we can log whether it was a successful exit or not
# code taken from: https://stackoverflow.com/a/9741784
class ExitHooks(object):
    def __init__(self):
        self.exit_code = None
        self.exception = None

    def hook(self):
        self._orig_exit = sys.exit
        sys.exit = self.exit
        self._orig_except = sys.excepthook
        sys.excepthook = self.exc_handler

    def exit(self, code=0):
        self.exit_code = str(code)
        self._orig_exit(code)

    def exc_handler(self, exc_type, exc, tb):
        self.exception = exc_type.__name__ + ": " + str(exc)
        self._orig_except(exc_type, exc, tb)

hooks = ExitHooks()
hooks.hook()


def logCmdLine( command ):

    global printLog, logFile

    logFile = "aligator.log"

    for idx,arg in enumerate(command):
        if re.search(r"(\s|\*)", arg):
            command[idx] = "'"+arg+"'"

    p = subprocess.Popen(['git', '--git-dir', os.path.dirname(command[0])+"/../.git",
    					  '--work-tree', os.path.dirname(command[0])+"/../",
                          'describe', '--always','--dirty','--tags'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    VERSION = p.communicate()[0].decode().strip()

    logStatement = "\n%s -- ALIGaToR %s run with command:\n\t%s\n" % (time.strftime("%c"), VERSION, " ".join(command))
    print(logStatement, file=sys.stderr)

    try:
        with open(logFile, "a") as handle:
            handle.write( logStatement )

        printLog = True

    except:
        print("Directory appears to be read-only; command line and output will not be saved", file=sys.stderr)


def logExit():

    global printLog, logFile
    if printLog:
        with open(logFile, "a") as handle:
            if hooks.exit_code is not None:
                formatted = re.sub( "\n", "\n\t", hooks.exit_code.strip(" \t\r\n") ) #remove white space and new lines on both ends; indent if multiple lines
                handle.write( "%s -- Program exited with error:\n\t%s\n" % (time.strftime("%c"),formatted) )
            elif hooks.exception is not None:
                formatted = re.sub( "\n", "\n\t", hooks.exception.strip(" \t\r\n") ) #remove white space and new lines on both ends; indent if multiple lines
                handle.write( "%s -- Exception:\n\t%s\n" % (time.strftime("%c"),formatted) )
            else:
                handle.write( "%s -- Program finished successfully\n" % time.strftime("%c") )

atexit.register(logExit)


########## UTILITIES ############

#make a codon table that can handle gaps
#start by getting the standard codon table
table = CodonTable.standard_dna_table.forward_table
#add gaps
for c1 in ["A", "C", "G", "T", "N"]:
	table["%s--"%c1] = "X"
	table["-%s-"%c1] = "X"
	table["--%s"%c1] = "X"
	for c2 in ["A", "C", "G", "T", "N"]:
		table["%s%s-"%(c1,c2)] = "X"
		table["-%s%s"%(c1,c2)] = "X"
		table["%s-%s"%(c1,c2)] = "X"
table["---"]="-"
#now register is and export
CodonTable.register_ncbi_table(name='gapped',alt_name="CAS0",id=99,table=table, stop_codons=['TAA', 'TAG', 'TGA', ], start_codons=['TTG', 'CTG', 'ATG', ] )
GAPPED_CODON_TABLE=CodonTable.ambiguous_dna_by_name["gapped"]



def load_seqs_in_dict(f, ids):
	result = dict()
	for entry in SeqIO.parse(open(f, "r"), "fasta"):
		if entry.id in ids:
			result[entry.id] = entry

	return result



def quickAlign( refseq, testseq, gapopen=None ):

	refseq = str( refseq.seq )
	refseq	= re.sub( "-", "", refseq )
	testseq = str( testseq.seq )
	testseq	= re.sub( "-", "", testseq )

	handle = StringIO()
	handle.write( ">ref\n%s\n>test\n%s\n"%(refseq,testseq) )
	data = handle.getvalue()

	muscle_cline = [ "muscle", '-quiet' ]
	if gapopen is not None: muscle_cline += [ '-gapopen', gapopen ]

	muscle_result = subprocess.run( muscle_cline, capture_output=True, text=True, input=data )

	aligned = dict()
	for p in SeqIO.parse(StringIO(muscle_result.stdout), "fasta"):
		aligned[ p.id ] = str(p.seq)
	return aligned
