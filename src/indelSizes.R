library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)

VC_dir <- args[1]
#VC_dir <- "/Users/mburger/Data/crispr/AFGHR/VariantCalls"
out_dir <- args[2]
#out_dir <- "/Users/mburger/Data/crispr/AFGHR/indel_quant/plots_trial"

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

#Plot histogram of all sizes detected
samples <- list.files(VC_dir)
widths <- list()
for (s in samples){
  genes <- list.files(file.path(VC_dir,s))
  for (g in genes){
    indel_events <- read.csv(file.path(VC_dir,s,g,"indel_events.csv"))
    indel_events <- indel_events[indel_events$onTarget == 1,]
    mult_vec <- rep(1,nrow(indel_events))
    mult_vec[indel_events[,"indel_type"] == "D"] <- -1
    indel_events$width <- indel_events$width * mult_vec
    df <- indel_events[,c("gene","width")]
    widths[[length(widths) + 1]] <- df
  }
}

df <- do.call(rbind,widths)
write.table(df,file=file.path(out_dir,"Indels_Observed.txt"),sep="\t",row.names=F,col.names=T)
pdf(file.path(out_dir,"Indels_Observed.pdf"),width=8,height=6)
#ggplot(df, aes(width, color = gene, fill=gene)) + 
print(ggplot(df, aes(width, color = gene)) + 
        geom_density(alpha = 0.3) + 
        scale_colour_Publication() +
        theme_bw() + 
        xlab("indel") + 
        scale_x_continuous(limits = c(-50, 50)) + 
        ggtitle("Observed Indel Widths") + 
        theme(panel.grid.major = element_blank(),panel.grid.major.x = element_blank()) + 
        theme(panel.grid.minor = element_blank()))
dev.off()
