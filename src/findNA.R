library("ggplot2")
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
# input_dir <- "/Users/mburger/Data/crispr/AGV28_QC/Power/VariantCalls"
out_dir <- args[2]
# out_dir <- "/Users/mburger/Data/crispr/AGV28_QC/Power"
max_width <- 150


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

#Determine which reads are NA
samples <- list.files(input_dir)
for (s in samples){
  sample_input_dir <- file.path(input_dir,s)
  genes <- list.files(sample_input_dir)
  df_list <- lapply(genes,function(g){
    WT_info <- read.csv(file.path(sample_input_dir,g,"all_reads_WT_info.csv"), stringsAsFactors=FALSE)
    WT_info <- WT_info[order(WT_info[,4],decreasing=T),]
    WT_info <- WT_info[!duplicated(WT_info[,1]),]
    if (nrow(WT_info) > 0){
      return(WT_info[is.na(WT_info[,"WT"]),])
    } else {
      return (list())
    }
  })
  df <- do.call(rbind,df_list)
  write.csv(df[,1],file=file.path(out_dir,paste0(s,"_qnames.csv")))
}




