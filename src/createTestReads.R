library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)

#WT_read_dir <- "/Users/mburger/Data/crispr/AGV28_QC/Power"
#indel_range <- 1:150
#out_fastq <- "/Users/mburger/Data/crispr/AGV28_QC/Power/tests.fastq"
#cutsite_annot  <- "~/Documents/crispr_dna/Annotations/cutsites.csv"
#gene_annot <- "/Users/mburger/Data/crispr/AGV28/AGV28_indel_genes.csv"
#read_length <- 300

WT_read_dir <- args[1]
indel_range <- 1:150
out_fastq <- file.path(WT_read_dir,"tests.fastq")
cutsite_annot  <- args[2]
read_length <- 300
gene_annot <- args[3]

indel_genes <- read.csv(gene_annot,header=F,stringsAsFactors = F)
rownames(indel_genes) <- indel_genes[,1]
#f1 <- FastqSampler("/Users/mburger/Data/crispr/AFVCC/1_AFVCC.1.1.fastq", n=50000)
#sample_reads <- yield(f1)
#writeFastq(sample_reads, fastq_file, "a",compress=FALSE)

#import reads for each gene and take 100 samples of each (if possible)
flist <- list.files(WT_read_dir,pattern=".txt")
gene_samples <- lapply(flist,function(f){
  tmp <- read.delim(file.path(WT_read_dir,f),stringsAsFactors=FALSE)
  if (nrow(tmp) > 100){
    tmp <- tmp[sample.int(nrow(tmp),100),]
  }
  return(tmp)
})

if (length(gene_samples) > 1){
  all_samples <- do.call(rbind,gene_samples)
} else {
  all_samples <- gene_samples[[1]]
}

#plot the number of samples for each gene
counts <- table(all_samples$SYMBOL)
counts <- sort(counts,decreasing=T)
cols <- rep("green",length(counts))
cols[counts < 100] <- "red"
pdf(file=file.path(WT_read_dir,"sample_sizes.pdf"),width=6,height=4)
barplot(counts,las=2,col=cols,ylab="WT Sample Reads")
dev.off()

#drop genes that don't have 100 WT samples
all_samples <- all_samples[all_samples$SYMBOL %in% names(counts)[counts == 100],,drop=F]

input_cutsite <- read.csv(cutsite_annot, header=FALSE, stringsAsFactors=FALSE)
rownames(input_cutsite) <- input_cutsite[,1]

#write WT reads
for (i in 1:nrow(all_samples)){
    
    seq <- all_samples[i,"SEQ"]
    qual <- all_samples[i,"QUAL"]
    pos_inc <- all_samples[i,"POS_INC"]
    qname <- paste0("@",all_samples[i,"QNAME"],":",0,":W:",NA)
    write(qname,file=out_fastq,append=T)
    write(seq,file=out_fastq,append=T)
    write("+",file=out_fastq,append=T)
    write(qual,file=out_fastq,append=T)
    
}

#create the deletion reads
for (i in 1:nrow(all_samples)){
  for (k in indel_range){
    
    seq <- all_samples[i,"SEQ"]
    qual <- all_samples[i,"QUAL"]
    pos_inc <- all_samples[i,"POS_INC"]
    cut_pos <- input_cutsite[all_samples[i,"SYMBOL"],3]
    
    if (k %% 2 == 0){
      side_drop <- k / 2
      str_cut_pos <- cut_pos - pos_inc
      str_cut_start <- str_cut_pos - side_drop
      str_cut_end <- str_cut_pos + side_drop
      seq <- paste0(substr(seq,1,str_cut_start),substr(seq,str_cut_end+1,nchar(seq)))
      qual <- paste0(substr(qual,1,str_cut_start),substr(qual,str_cut_end+1,nchar(qual)))
    } else {
      side_drop <- (k - 1) / 2
      str_cut_pos <- cut_pos - pos_inc
      str_cut_start <- str_cut_pos - side_drop
      str_cut_end <- str_cut_pos + side_drop
      seq <- paste0(substr(seq,1,str_cut_start),substr(seq,str_cut_end+2,nchar(seq)))
      qual <- paste0(substr(qual,1,str_cut_start),substr(qual,str_cut_end+2,nchar(qual)))
    }
    
    if (all_samples[i,"READ_DIR"] == "+"){
      seq <- paste0(seq,paste(rep("A",k),collapse=""))
      qual <- paste0(qual,paste(rep("#",k),collapse=""))
    } else {
      seq <- paste0(paste(rep("A",k),collapse=""),seq)
      qual <- paste0(paste(rep("#",k),collapse=""),qual)
    }
    
    indel_pos <- cut_pos - side_drop
    qname <- paste0("@",all_samples[i,"QNAME"],":",k,":D:",indel_pos)
    write(qname,file=out_fastq,append=T)
    write(seq,file=out_fastq,append=T)
    write("+",file=out_fastq,append=T)
    write(qual,file=out_fastq,append=T)
  }
  
}

#create the insertion reads
for (i in 1:nrow(all_samples)){
  for (k in indel_range){
    
    symbol <- all_samples[i,"SYMBOL"]
    
    seq <- all_samples[i,"SEQ"]
    qual <- all_samples[i,"QUAL"]
    pos_inc <- all_samples[i,"POS_INC"]
    cut_pos <- input_cutsite[symbol,3]
    amplicon_length <- indel_genes[symbol,4] - indel_genes[symbol,3] + 1
    
    str_cut_pos <- cut_pos - pos_inc
    splt <- unlist(strsplit(seq,""))
    
    ins <- paste(splt[sample.int(nchar(seq),size=k)],collapse="")
    seq <- paste0(substr(seq,1,str_cut_pos),ins,substr(seq,str_cut_pos+1,nchar(seq)))
    seq <- substr(seq,1,read_length)
    
    qual_ins <- paste(rep('5',k),collapse="")
    qual <- paste0(substr(qual,1,amplicon_length),qual_ins,substr(qual,amplicon_length+1,nchar(qual)))
    qual <- substr(qual,1,read_length)
    
    qname <- paste0("@",all_samples[i,"QNAME"],":",k,":I:",cut_pos)
    write(qname,file=out_fastq,append=T)
    write(seq,file=out_fastq,append=T)
    write("+",file=out_fastq,append=T)
    write(qual,file=out_fastq,append=T)
  }
  
}

