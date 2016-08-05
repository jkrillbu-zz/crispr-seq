library("ShortRead")

#fqQC = qa(dirPath="/Users/mburger/Data/crispr/tmp/", pattern=".fastq$", type="fastq")
#report(fqQC, type="html", dest="/Users/mburger/Data/crispr/tmp/fastqQAreport")

args <- commandArgs(trailingOnly = TRUE)

output_folder <- args[1]
if (!file.exists(output_folder)){
  dir.create(output_folder,recursive = T)
}
setwd(output_folder)

fl_barcodes <- args[2]
fl_reads <- args[3]
barcode_file <- args[4]

reads_folder <- "Reads"
dir.create(reads_folder)

barcodes <- read.csv(barcode_file, header=FALSE, stringsAsFactors=FALSE)
samples <- barcodes[,2]
names(samples) <- barcodes[,1]
samples <- samples[samples != ""]
names(samples) <- gsub("[^A-Za-z0-9_]","",names(samples))
sample_files <- paste0(names(samples),".fastq")
counts <- rep(0,length(samples)+1)

f_barcodes <- FastqStreamer(fl_barcodes, 50000)
f_reads <- FastqStreamer(fl_reads, 50000)

pdict <- PDict(samples)

progress <- 0

repeat {
  
  fq_barcodes <- yield(f_barcodes)
  if (length(fq_barcodes) == 0)
    break
  fq_reads <- yield(f_reads)
  
  fq_barcodes <- trimTails(object=fq_barcodes, k=2, a="?")
  trimmed <- width(fq_barcodes) < 8
  fq_barcodes <- fq_barcodes[!trimmed]
  fq_reads <- fq_reads[!trimmed]
  
  
  barcodes <- sread(fq_barcodes)
  mappings <- vwhichPDict(pdict,barcodes,max.mismatch=0)
  
  for (i in 1:length(mappings)){
    if (length(mappings[[i]]) != 1) {
      mappings[[i]] <- 0
    }
  }
  
  mappings <- unlist(mappings)
  
  for (i in 1:length(sample_files)){
    
    filter <- mappings == i
    counts[i] <- counts[i] + sum(filter)
    writeFastq(fq_reads[filter], file.path(reads_folder,sample_files[i]), "a",compress=FALSE)
    
  }
  
  filter <- mappings == 0
  counts[length(counts)] <- counts[length(counts)] + sum(filter)
  writeFastq(fq_reads[filter], "unknown_reads.fastq", "a",compress=FALSE)
  writeFastq(fq_barcodes[filter], "unknown_barcodes.fastq", "a",compress=FALSE)
  
  progress <- progress + 50000
  print(progress)
}

sample_ids <- c(names(samples),"unknown")
barcode_list <- c(samples,"********")
barcode_mapping <- data.frame(cbind(sample_ids,barcode_list,counts))
write.table(barcode_mapping,file="barcode_mapping.csv",sep=",",row.names=F,col.names=T)
