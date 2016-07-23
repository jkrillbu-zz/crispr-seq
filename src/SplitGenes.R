#import the genes that were targeted for sequencing 

args <- commandArgs(trailingOnly = TRUE)
indel_genes_file <- args[1]
reads_folder <- args[2]

gene_folder <- file.path(reads_folder,"gene_level")
if (!file.exists(gene_folder)){
  dir.create(gene_folder,recursive = T)
}

indel_genes <- read.csv(indel_genes_file, header=FALSE, stringsAsFactors=FALSE)
indel_genes <- indel_genes[indel_genes[,5],]

filelist <- list.files(reads_folder,pattern=".sorted.bam$")

samtools <- "samtools view -b"

for (f in filelist){
  sample <- gsub(".sorted.bam$","",f)
  for (i in 1:nrow(indel_genes)){
    filename <- paste0(sample,"_",indel_genes[i,1],".bam")
    fp <- file.path(gene_folder,filename)
    location <- paste0(indel_genes[i,2],":",indel_genes[i,3],"-",indel_genes[i,4])
    cmd <- paste(samtools,file.path(reads_folder,f),location,">",fp,sep=" ")
    system(cmd)
  }
}