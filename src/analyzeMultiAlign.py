import re
import csv
import pysam
import sys, getopt
import os
import glob
from subprocess import call

#1 indexed
def matchEndRead(cig):

	length = 0

	for (e,l) in cig:
		if e == 0:
			return length + l
		else:
			length += l

	return None

#1 indexed
def matchStartRead(cig):

	length = 0

	for (e,l) in cig:
		if e == 0:
			return length + 1
		else:
			length += l

	return None

def matchEnd(pos,cig):

	length = 0

	for (e,l) in cig:
		if e == 0 or e == 1 or e == 2:
			length += l

	return (pos + length - 1)

def alignmentScore(tags):

	for (e,l) in tags:
		if e == 'AS':
			return l

	return 0

def indelPosition(left_read,right_read):
	left_pos = matchEnd(left_read.pos,left_read.cigar) + 1
	right_pos = right_read.pos - 1

	return ((left_pos,right_pos))

def indelPositionRead(left_read,right_read):
	left_pos = matchEndRead(left_read.cigar) + 1
	right_pos = matchStartRead(right_read.cigar) - 1

	return((left_pos,right_pos))


def findPair(read_list,target_interval,gene_interval):

	left_read = None
	left_score = 0
	right_read = None
	right_score = 0

	for read in read_list:
		AS = alignmentScore(read.tags)
		if (read.pos >= gene_interval[1] - 5) and matchEnd(read.pos,read.cigar) <= target_interval[1] and AS > left_score:
			left_read = read
			left_score = AS
		elif read.pos >= target_interval[1] and (matchEnd(read.pos,read.cigar) <= gene_interval[2] + 5) and AS > right_score:
			right_read = read
			right_score = AS

	return((left_read,right_read))

#Writes the width of the gap in the read for insertions and the width of the gap of the reference for deletions
#Complex combinations of insertion and deletion both get written
def analyzeMult(sample,gene,read_list,target_interval,gene_interval,WT_writer,fd_writer):

	if len(read_list) > 1:
		left_read,right_read = findPair(read_list,target_interval,gene_interval)
		if (left_read is not None) and (right_read is not None):
			left_pos,right_pos = indelPosition(left_read,right_read)
			total_length = matchEnd(right_read.pos,right_read.cigar) - left_read.pos + 1
			row = [left_read.qname,str(total_length),'0','4']
			WT_writer.writerow(row)
			width = right_pos - left_pos + 1
			left_pos_read,right_pos_read = indelPositionRead(left_read,right_read)
			read_width = right_pos_read - left_pos_read + 1
			
			if read_width < width:
				row = [sample,gene,left_read.qname,str(total_length),'D',str(left_pos+1),str(right_pos+1),str(width),'1']
				fd_writer.writerow(row)
			elif read_width > width:
				row = [sample,gene,left_read.qname,str(total_length),'I',str(left_pos+1),str(right_pos+1),str(read_width),'1']
				fd_writer.writerow(row)
			else:
				row = [sample,gene,left_read.qname,str(total_length),'D',str(left_pos+1),str(right_pos+1),str(width),'1']
				fd_writer.writerow(row)
				row = [sample,gene,left_read.qname,str(total_length),'I',str(left_pos+1),str(right_pos+1),str(read_width),'1']
				fd_writer.writerow(row)

		#else 
			#write reads to problem bam
	#write read to problem bam
	return(1)


def iterateReads(sample,gene,bam_file,target_interval,gene_interval,files):

	WT_fd = open(files['WT'],'a')
	WT_writer = csv.writer(WT_fd)

	fd = open(files['indel'],'a')
	fd_writer = csv.writer(fd)

    #Sort the NA.bam which is a file full of the reads that didn't span the cut interval
	#pysam.sort('-n', bam_file, os.path.join(os.path.dirname(bam_file),'NA.sorted'), catch_stdout=False)
	call(["samtools", "sort", "-n", bam_file, os.path.join(os.path.dirname(bam_file),"NA.sorted")])
	bam_file = os.path.join(os.path.dirname(bam_file),'NA.sorted.bam')

	bamFP = pysam.Samfile(bam_file, "rb")
	#iter1 = bamFP.fetch(until_eof=True)
	#read = iter1.next()
	read_list = []
	qname = ''
	
	for read in bamFP.fetch(until_eof=True):
		if len(read_list) == 0 or read.qname == qname:
			read_list.append(read)
			qname = read.qname
		else:
			a = analyzeMult(sample,gene,read_list,target_interval,gene_interval,WT_writer,fd_writer)
			read_list = [read]
			qname = read.qname

	#analyze last read_list

	bamFP.close()
	fd.close()
	WT_fd.close()
	

	return(1)


def main(argv):

	inputfile = ''
	genes = ''
	cutsites = ''
	try:
		opts, args = getopt.getopt(argv,"h:i:c:g:")
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -c <cutsites> -g <genes>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -c <cutsites> -g <genes>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-c", "--cutsites"):
			cutsites = arg
		elif opt in ("-g", "--genes"):
			genes = arg
	

	gene_dict = {}

	with open(genes, 'rU') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in spamreader:
			gene_dict[row[0]] = [row[1],int(row[2])-1,int(row[3])-1]

	target_dict = {}
	
	with open(cutsites, 'rU') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in spamreader:
			target_dict[row[0]] = [int(row[2])-1+20,int(row[3])-1-25]


	flist = glob.glob(os.path.join(inputfile,"*"))

	for f in flist:

		sample = os.path.basename(f)

		#NA_outdir = os.path.dirname(f)
		#NA_outfile = os.path.join(NA_outdir,sample + '_frag.bam')

		#samfile = pysam.AlignmentFile(f, "rb")
		#problemBAM = pysam.AlignmentFile(NA_outfile, "wb", template=samfile)
		#samfile.close()

		for gene in gene_dict:
			outdir = os.path.join(f,gene)

			files = {}
			files['WT'] = os.path.join(outdir,'all_reads_WT_info.csv')
			files['indel'] = os.path.join(outdir,'indel_events.csv')
			
			bam_file = os.path.join(outdir,'NA.bam')

			z = iterateReads(sample,gene,bam_file,target_dict[gene],gene_dict[gene],files)

		#problemBAM.close()

	return(1)


if __name__ == "__main__":
   main(sys.argv[1:])
