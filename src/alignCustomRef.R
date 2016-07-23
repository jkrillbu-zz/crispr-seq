args <- commandArgs(trailingOnly = TRUE)

reads_dir <- args[1]
custom_ref <- args[2]
known_indel_samples_annot <- args[3]
out_file <- args[4]


known_indel_samples <- read.csv(known_indel_samples_annot, stringsAsFactors=FALSE)
known_indel_samples <- known_indel_samples[known_indel_samples[,1] != "",]
known_indel_samples[,1] <- gsub("[^A-Za-z0-9_]","",known_indel_samples[,1])
rownames(known_indel_samples) <- known_indel_samples[,1]
known_indel_samples <- known_indel_samples[,c(-1),drop=F]

r <- nrow(known_indel_samples)
c <- ncol(known_indel_samples)
output <- matrix(data=rep("",r*c),nrow=r,ncol=c)
rownames(output) <- rownames(known_indel_samples)
colnames(output) <- colnames(known_indel_samples)

for (sample in rownames(known_indel_samples)){
  cmd <- paste0("/Users/mburger/Documents/bwa/bwa mem -M -t 6 ",custom_ref," ",file.path(reads_dir,paste0(sample,".fastq"))," > ",file.path(reads_dir,paste0(sample,".custom.sam")))
  system(cmd)
  
  cmd <- paste0("/usr/local/bin/samtools view -q 4 ",file.path(reads_dir,paste0(sample,".custom.sam"))," > ",file.path(reads_dir,paste0(sample,"_mapped.custom.sam")))
  system(cmd)
  system(paste0("rm ",file.path(reads_dir,paste0(sample,".custom.sam"))))
  
  no_col <- max(count.fields(file.path(reads_dir,paste0(sample,"_mapped.custom.sam")), sep = "\t"))
  if (no_col != -Inf){
    CD19 <- read.table(file.path(reads_dir,paste0(sample,"_mapped.custom.sam")), header = FALSE,sep="\t",fill=TRUE,colClasses = c("character","integer","character","integer","integer","character",rep("NULL", no_col - 6)),col.names=c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR",as.character(7:no_col)))
  }
  genes <- colnames(known_indel_samples)[known_indel_samples[sample,]]
  for (g in genes){
    if (no_col != -Inf){
      WT <- paste0(g,"_WT")
      MUT <- paste0(g,"_MUT")
      
      RNAME <- CD19[CD19[,"RNAME"] %in% c(WT,MUT),"RNAME"]
      total_reads <- length(RNAME)
      RNAME <- RNAME[RNAME %in% c(MUT)]
      mut_reads <- length(RNAME)
    } else {
      mut_reads <- 0
      total_reads <- 0
    }
    
    if (total_reads == 0 || mut_reads == 0){
      frac <- 0
    } else {
      frac <- mut_reads / total_reads
    }
    output[sample,g] <- paste0(frac,";",total_reads)
  }
  
}

write.csv(output,file=out_file)