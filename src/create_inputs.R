
crispr_dir <- '/Users/mburger/Data/crispr'
project_name <- 'AU6R0'

gRNAs <- read.csv("~/Data/crispr/AU6R0/gRNAs.csv", stringsAsFactors=FALSE)
gRNAs <- gRNAs[gRNAs[,1] != "",]

#indel_genes
amplicon_split <- strsplit(gRNAs$PCR.amplicon,'[:-]')
indel_genes <- data.frame(Gene=gRNAs$Gene,chr=sapply(amplicon_split,function(x){x[1]}),
                          start=sapply(amplicon_split,function(x){x[2]}),
                          end=sapply(amplicon_split,function(x){x[3]}),
                          include=rep(TRUE,nrow(gRNAs)),strand=gRNAs$Strand)
  
write.table(indel_genes,file=file.path(crispr_dir,project_name,paste0(project_name,'_indel_genes.csv')),sep=",",col.names = F, row.names=F,quote = F)

#inel_cut_site
cut_split <- strsplit(gRNAs$Cut.site,':')
indel_cut_site <- indel_genes[,1:2]
indel_cut_site <- cbind(indel_cut_site,sapply(cut_split,function(x){x[2]}))
write.table(indel_cut_site,file=file.path(crispr_dir,project_name,paste0(project_name,'_indel_cut_site.csv')),sep=",",col.names = F, row.names=F,quote = F)

#indel_cut_interval
cut_sites <- as.numeric(as.character(indel_cut_site[,3]))
indel_cut_interval <- indel_genes[,1:2]
indel_cut_interval <- cbind(indel_cut_interval,cut_sites - 25,cut_sites + 25)
write.table(indel_cut_interval,file=file.path(crispr_dir,project_name,paste0(project_name,'_indel_cut_interval.csv')),sep=",",col.names = F, row.names=F,quote = F)
  
#negative_controls
negative_controls <- read.table("~/Data/crispr/AU6R0/control.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
if (ncol(negative_controls) == 1){
  negative_controls <- matrix(rep(1,nrow(negative_controls)*nrow(gRNAs)),nrow=nrow(negative_controls),ncol=nrow(gRNAs),dimnames=list(negative_controls[,1],gRNAs$Gene))
}
write.csv(negative_controls,file=file.path(crispr_dir,project_name,paste0(project_name,'_negative_controls.csv')),quote=F)
  