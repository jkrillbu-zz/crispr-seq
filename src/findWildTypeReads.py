import re
import csv
import pysam
import sys, getopt
import os
import glob

# @param cig: cigar string object from read
# return: True if any parts of the cigar string are D (2) or I (1)
def hasIndel(cig):

	for (e,l) in cig:
		if e == 1 or e == 2 or e == 5:
			return(True)

	return(False)

# @param cig: cigar string object from read
# @param a_len: actual length of PCR amplicon
# return: True if M (0), match and mismatch, part of cigar string is length of amplicon
def isMatch(cig,a_len):

	for (e,l) in cig:
		if e == 0 and l >= a_len:
			return(True)

	return(False)

# @param ts: the list of tags
# return: the number of mismatches
def countMismatches(ts):

	for (e,l) in ts:
		if e == 'NM':
			return(l)

	return(0)

# @param cig: cigar string
# @param pos: start position of non-soft clipped read
# return: start position of read including soft clipped prefix
def findSeqStart(cig,pos):

	for (e,l) in cig:
		if e == 4:
			pos = pos - l
		elif e == 0:
			break

	return(pos)

def runGene(inputfile,outputfile,gene,gene_interval):

	print([inputfile,outputfile,gene,gene_interval])

	WT_fd = open(outputfile,'a')
	WT_writer = csv.writer(WT_fd)
	
	flist = glob.glob(os.path.join(inputfile,"*.sorted.bam"))

	for f in flist:
		count = 0
		hasNext = True
		bamFP = pysam.AlignmentFile(f, "rb")
		iter1 = bamFP.fetch(gene_interval[0], gene_interval[1], gene_interval[2])

		while (hasNext and (count < 100)):
			try: read = iter1.next()
			except StopIteration: hasNext = False
			else: hasNext = True
			if hasNext:
				start_match = read.pos <= (gene_interval[1] - 1)
				end_match = read.reference_end >= gene_interval[2]
				#seq_length = read.query_alignment_length == 300
				seq_length = True
				no_indel = not(hasIndel(read.cigar))
				good_quality = read.mapq == 60
				match_length = isMatch(read.cigar,(gene_interval[2]-gene_interval[1]+1))
				if (start_match and end_match and seq_length and no_indel and good_quality and match_length):
					count = count + 1
					refchar = read.rname + 1
					if refchar == 23:
						refchar = 'X'
					elif refchar == 24:
						refchar = 'Y'
					read_direction = "+"
					if read.is_reverse:
						read_direction = "-" 
					row = [gene + "\t" + read.qname + "\t" + str(refchar) + "\t" + str(findSeqStart(read.cigar,read.pos)+1) + "\t" + str(read.pos+1) + "\t" + read.cigarstring + "\t" + read.seq + "\t" + read.qual + "\t" + str(countMismatches(read.tags)) + "\t" + read_direction]
					WT_writer.writerow(row)

	
		bamFP.close()
	WT_fd.close()

	return(count)

def main(argv):

	inputfile = ''
	outputdir = ''
	genes = ''

	try:
		opts, args = getopt.getopt(argv,"h:i:o:g:")
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -o <outputdir> -g <genes>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputdir> -g <genes>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputdir = arg
		elif opt in ("-g", "--genes"):
			genes = arg

	gene_dict = {}

	with open(genes, 'rU') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in spamreader:
			gene_dict[row[0]] = [row[1],int(row[2]),int(row[3])]


	for g in gene_dict:
		outputfile = os.path.join(outputdir,"sample_" + g + ".txt")
		WT_fd = open(outputfile,'a')
		WT_writer = csv.writer(WT_fd)
		row = ["SYMBOL\tQNAME\tRNAME\tPOS_INC\tPOS\tCIGAR\tSEQ\tQUAL\tNM\tREAD_DIR"]
		WT_writer.writerow(row)
		WT_fd.close()
		res = runGene(inputfile,outputfile,g,gene_dict[g])

	return(1)


if __name__ == "__main__":
   main(sys.argv[1:])
