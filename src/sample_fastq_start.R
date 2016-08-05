library("ShortRead")

#fqQC = qa(dirPath="/Users/mburger/Data/crispr/tmp/", pattern=".fastq$", type="fastq")
#report(fqQC, type="html", dest="/Users/mburger/Data/crispr/tmp/fastqQAreport")

args <- commandArgs(trailingOnly = TRUE)

#output_folder <- args[1]
#if (!file.exists(output_folder)){
#  dir.create(output_folder,recursive = T)
#}
#setwd(output_folder)

fl_barcodes <- args[1]
fl_reads <- args[2]
max_reads <- as.numeric(args[3])
run_id <- args[4]

f_barcodes <- FastqStreamer(fl_barcodes, 50000)
f_reads <- FastqStreamer(fl_reads, 50000)

progress <- 0

while (progress < max_reads) {
  
  fq_barcodes <- yield(f_barcodes)
  if (length(fq_barcodes) == 0)
    break
  fq_reads <- yield(f_reads)
  
  writeFastq(fq_barcodes, paste0(run_id,"_sample_barcodes.fastq"), "a",compress=FALSE)
  writeFastq(fq_reads, paste0(run_id,"_sample_reads.fastq"), "a",compress=FALSE)

  progress <- progress + 50000
  print(progress)
}
