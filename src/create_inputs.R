
args <- commandArgs(trailingOnly = TRUE)
gRNA_file <- args[1]
#gRNA_file <- "/Users/mburger/dockeryard/test_data/AU6R0_sample_gRNAs.csv"
neg_control_file <- args[2]
#neg_control_file <- "~/Data/crispr/AU6R0/control.csv"

gRNAs <- read.csv(gRNA_file, stringsAsFactors=FALSE)
gRNAs <- gRNAs[gRNAs[,1] != "",]

#indel_genes
amplicon_split <- strsplit(gRNAs$amplicon,'[:-]')
indel_genes <- data.frame(Gene=gRNAs$gene,chr=sapply(amplicon_split,function(x){x[1]}),
                          start=sapply(amplicon_split,function(x){x[2]}),
                          end=sapply(amplicon_split,function(x){x[3]}),
                          include=rep(TRUE,nrow(gRNAs)),strand=gRNAs$strand)
  
write.table(indel_genes,file='indel_genes.csv',sep=",",col.names = F, row.names=F,quote = F)

#inel_cut_site
cut_split <- strsplit(gRNAs$cut,':')
indel_cut_site <- indel_genes[,1:2]
indel_cut_site <- cbind(indel_cut_site,sapply(cut_split,function(x){x[2]}))
write.table(indel_cut_site,file='indel_cut_site.csv',sep=",",col.names = F, row.names=F,quote = F)

#indel_cut_interval
cut_sites <- as.numeric(as.character(indel_cut_site[,3]))
indel_cut_interval <- indel_genes[,1:2]
indel_cut_interval <- cbind(indel_cut_interval,cut_sites - 25,cut_sites + 25)
write.table(indel_cut_interval,file='indel_cut_interval.csv',sep=",",col.names = F, row.names=F,quote = F)
  
#negative_controls
negative_controls <- read.csv(neg_control_file, row.names=1, stringsAsFactors=FALSE)
if (ncol(negative_controls) < 1){
  negative_controls <- matrix(rep(1,nrow(negative_controls)*nrow(gRNAs)),nrow=nrow(negative_controls),ncol=nrow(gRNAs),dimnames=list(rownames(negative_controls),gRNAs$gene))
}
write.csv(negative_controls,file='negative_controls.csv',quote=F)
  