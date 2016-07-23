import re
import csv
import pysam
import sys, getopt
import os
import glob

def parse(cig,start):

	results = []

	for num1, i_or_d, num2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', cig):
		results.append((num1,i_or_d,str(start)))
		if num1:
			start += int(num1)
		if num2:
			start += int(num2)

	return(results)

def parse2(cig,start):

	results = []
	ref_pos = start

	for (e,l) in cig:
		if e == 1 or e == 2:
			#add to results
			if e == 1:
				results.append((l,'I',ref_pos))
				ref_pos += 1
			else:
				results.append((l,'D',ref_pos))
				ref_pos += l
		elif e == 0:
			ref_pos += l

	return(results)

def alignmentScore(tags):

	for (e,l) in tags:
		if e == 'AS':
			return l

	return 0

def alignedReadLength(cig):

	length = 0

	for (e,l) in cig:
		if e == 0 or e == 1 or e == 2:
			length += l

	return length

def detectOverlap(target_interval,indel_interval):

	if indel_interval[0] <= target_interval[1] and indel_interval[1] >= target_interval[0]:
		return True
	else:
		return False

def isWT(result_list,target_interval):

	WT = True

	if len(result_list) == 0:
		return(WT)

	for (width,event,start) in result_list:
		end = start + width - 1
		if detectOverlap(target_interval,[start,end]):
			WT = False
			break

	return(WT)

def printIndels(result_list,sample,gene,read_id,read_length,target_interval,fd_writer):

	for (width,event,start) in result_list:
		end = start + width - 1
		if detectOverlap(target_interval,[start,end]):
			on_target = 1
		else:
			on_target = 0
		row = [sample,gene,read_id,read_length,event,str(start+1),str(end+1),width,on_target]
		fd_writer.writerow(row)

	return(1)

def printWT(result_list,read_id,start,length,files,target_interval,tags,WT_writer,problemBAM,read):

	if isWT(result_list,target_interval):
		WT = '1'
	else:
		WT = '0'

	end = start + length - 1

	if not(start <= target_interval[0] and end >= target_interval[1]):
		WT = 'NA'
		problemBAM.write(read)

	row = [read_id,length,WT,alignmentScore(tags)]
	WT_writer.writerow(row)

	return(1)

def iterateReads(bam_file,sample,gene,gene_interval,target_interval,files):

	WT_fd = open(files['WT'],'a')
	WT_writer = csv.writer(WT_fd)

	fd = open(files['indel'],'a')
	fd_writer = csv.writer(fd)

	samfile = pysam.AlignmentFile(bam_file, "rb")
	problemBAM = pysam.AlignmentFile(files['NA'], "wb", template=samfile)
	samfile.close()

	indel_header = ["sample","gene","read_id","read_length","indel_type","start","end","width","onTarget"]
	fd_writer.writerow(indel_header)

	wt_header = ["read_id","read_length","WT","AS"]
	WT_writer.writerow(wt_header)

	bamFP = pysam.AlignmentFile(bam_file, "rb")
	read_id = 1
	
	for read in bamFP.fetch(gene_interval[0], gene_interval[1], gene_interval[2]):
		#if not(read.is_secondary):
		start = read.pos
		length = alignedReadLength(read.cigar)
		r = parse2(read.cigar,start)
		x = printWT(r,read.qname,start,length,files,target_interval,read.tags,WT_writer,problemBAM,read)
		if len(r) > 0:
			y = printIndels(r,sample,gene,read.qname,length,target_interval,fd_writer)
		read_id += 1
			

	bamFP.close()
	fd.close()
	WT_fd.close()
	problemBAM.close()

	return(1)

#def runGene(sample,gene,basedir,inputfile,gene_interval,target_interval,problemBAM_writer):

#	files = {}
#	files['WT'] = os.path.join(basedir,'all_reads_WT_info.csv')
#	files['indel'] = os.path.join(basedir,'indel_events.csv')

#	z = iterateReads(inputfile,sample,gene,gene_interval,target_interval,files,problemBAM_writer)

#	return(1)

def main(argv):

	inputfile = ''
	outputfile = ''
	cutsites = ''
	try:
		opts, args = getopt.getopt(argv,"h:i:o:c:g:")
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -o <outputfile> -c <cutsites> -g <genes>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputfile> -c <cutsites> -g <genes>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
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
			target_dict[row[0]] = [int(row[2])-1,int(row[3])-1]


	flist = glob.glob(os.path.join(inputfile,"*.sorted.bam"))

	for f in flist:

		bam_file = os.path.basename(f)
		sample,t1,t2 = bam_file.split(".")

		#NA_outdir = os.path.dirname(f)
		#NA_outfile = os.path.join(NA_outdir,sample + '_NA.bam')

		#samfile = pysam.AlignmentFile(f, "rb")
		#problemBAM = pysam.AlignmentFile(NA_outfile, "wb", template=samfile)
		#samfile.close()

		for gene in gene_dict:
			outdir = os.path.join(outputfile, sample, gene)
			os.makedirs(outdir)
			files = {}
			files['WT'] = os.path.join(outdir,'all_reads_WT_info.csv')
			files['indel'] = os.path.join(outdir,'indel_events.csv')
			files['NA'] = os.path.join(outdir,'NA.bam')
			#e = runGene(sample,gene,outdir,f,gene_dict[gene],target_dict[gene],problemBAM)
			z = iterateReads(f,sample,gene,gene_dict[gene],target_dict[gene],files)

		#problemBAM.close()


	return(1)


if __name__ == "__main__":
   main(sys.argv[1:])
