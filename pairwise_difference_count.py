#!/usr/bin/env python

#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys
from optparse import OptionParser
from Bio.Align import AlignInfo
from Bio import AlignIO


def DoError(ErrorString):
        print "!!!Error:", ErrorString,"!!!"
        sys.exit()



def read_alignment(filename, quiet=False):
        if not quiet:
                print "Reading alignment file..."

        filetype=["phylip", "fasta", "clustal", "nexus", "emboss", "stockholm", "fasta-m10", "ig"]

        if filename.split(".")[-1].lower() in filetype:
                guesstype=filename.split(".")[-1].lower()
        elif filename.split(".")[-1].lower() in ["phy"]:
                guesstype="phylip"
        elif filename.split(".")[-1].lower() in ["fna", "dna", "aa", "aln", "fas"]:
                guesstype="fasta"
        elif filename.split(".")[-1].lower() in ["nxs", "nex", "nexus"]:
                guesstype="nexus"
        else:
                guesstype=""
        
        readok=False
        
        
        if guesstype!="":
                if not quiet:
                        print "Guessing file is in "+guesstype+" format"
                        print "Trying to open file "+filename+" as "+guesstype
                try:
                        alignmentObject = AlignIO.read(open(filename, "rU"), guesstype)
                except StandardError:
                        print "Cannot open alignment file as "+guesstype
                else:
                        readok=True
                        
                filetype.remove(guesstype)

        x=0

        while readok==False and x<len(filetype):
                if not quiet:
                        print "Trying to open file "+filename+" as "+filetype[x]

                try:
                        alignmentObject = AlignIO.read(open(filename), filetype[x])
                except StandardError:
                        print "Cannot open alignment file "+filename+" as "+filetype[x]
                else:
                        readok=True

                x=x+1

        if readok==False:
                raise SimonError("Failed to read alignment")
        else:
                if not quiet:
                        print "Alignment read successfully"
                return alignmentObject


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')


gap_and_missing=set(["-", "N", "?"])

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	#Read the alignment file
	
	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")

	seqnames=[]
	for sequence in alignment:
		if not sequence.name in seqnames:
			seqnames.append(sequence.name)
	
	print '\t'.join(["Taxon1", "Taxon2", "SNPs", "%ID", "Aligned bases"])
	for x, taxon in enumerate(seqnames):
		taxonseq=""
		for sequence in alignment:
			if sequence.name==taxon:
				taxonseq=str(sequence.seq).upper()
				break
		if taxonseq=="":
			print "Cannot find ", taxon
			continue
		for taxonb in seqnames[x+1:]:
						
			if taxon==taxonb:
				continue
			
			taxonbseq=""
			for sequence in alignment:
				if sequence.name==taxonb:
					taxonbseq=str(sequence.seq).upper()
					break
			
			if taxonbseq=="":
				print "Cannot find ", taxonb
				continue
			
			count=0
			ncount=0
			idcount=0			
			
			if len(taxonseq)==len(taxonbseq):
			
				for y, base in enumerate(taxonseq):
					
					if base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						idcount+=1
					
					if base!=taxonbseq[y] and base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						count+=1
						#print y+1, base, taxonbseq[y]
					elif base=="N":
						ncount+=1
						
			else:
				print taxon, "and", taxonb, "are not the same length"
#			pairwise_distances[taxonb][taxon]=count
			print '\t'.join(map(str, [taxon, taxonb, count, ((float(idcount)-float(count))/idcount)*100, idcount]))
