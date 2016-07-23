library("ggplot2")
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
#input_dir <- "/Users/mburger/Data/crispr/AFGHR/Power/VariantCalls/tests"
out_dir <- args[2]
#out_dir <- "/Users/mburger/Data/crispr/AFGHR/Power/VariantCalls"
max_width <- 150

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#026910","#BF9202")), ...)
  
}

#Plot the WT accuracy
genes <- list.files(input_dir)
df_list <- lapply(genes,function(g){
  WT_info <- read.csv(file.path(input_dir,g,"all_reads_WT_info.csv"), stringsAsFactors=FALSE)
  WT_info <- WT_info[order(WT_info[,4],decreasing=T),]
  WT_info <- WT_info[!is.na(WT_info[,"WT"]),]
  WT_info <- WT_info[!duplicated(WT_info[,1]),]
  if (nrow(WT_info) > 0){
    qname_split <- strsplit(WT_info[,1],":")
    q <- do.call(rbind,qname_split)
    WT_info <- cbind(q,WT_info[,2:3],stringsAsFactors=FALSE)
    
    mult_vec <- rep(1,nrow(WT_info))
    mult_vec[WT_info[,9] == "D"] <- -1
    WT_info[,8] <- as.numeric(WT_info[,8])
    WT_info[,8] <- WT_info[,8] * mult_vec
    WT_info[,9] <- as.numeric(WT_info[,9] == 'W')
    WT_info[is.na(WT_info[,"WT"]),"WT"] <- -1
    WT_info <- WT_info[WT_info[,"WT"] == WT_info[,9],]
    accuracy <- table(factor(as.character(WT_info[,8]),levels=-max_width:max_width))
    return(data.frame(indel=as.numeric(names(accuracy)),accuracy=as.numeric(accuracy),gene=g))
  } else {
    return (list())
  }
})
df <- do.call(rbind,df_list)
write.table(df,file=file.path(out_dir,"WT_Accuracy.txt"),sep="\t",row.names=F,col.names=T)
p <- ggplot(df, aes(x=indel, y=accuracy, group=gene,colour = gene))
pdf(file.path(out_dir,"WT_Accuracy.pdf"),width=8,height=6)
print(p + ggtitle("Test Set Accuracy") + geom_line(size=1) + scale_colour_Publication() + theme_bw() + scale_x_continuous(breaks = seq(-150, 150, by = 20)))
dev.off()

# #Plot histogram of all sizes detected
# VC_dir <- "/Users/mburger/Data/crispr/AFGHR/VariantCalls"
# samples <- list.files(VC_dir)
# widths <- list()
# for (s in samples){
#   genes <- list.files(file.path(VC_dir,s))
#   for (g in genes){
#     indel_events <- read.csv(file.path(VC_dir,s,g,"indel_events.csv"))
#     indel_events <- indel_events[indel_events$onTarget == 1,]
#     mult_vec <- rep(1,nrow(indel_events))
#     mult_vec[indel_events[,"indel_type"] == "D"] <- -1
#     indel_events$width <- indel_events$width * mult_vec
#     df <- indel_events[,c("gene","width")]
#     widths[[length(widths) + 1]] <- df
#   }
# }
# 
# df <- do.call(rbind,widths)
# pdf(file.path(out_dir,"Indels_Observed.pdf"),width=8,height=6)
# #ggplot(df, aes(width, color = gene, fill=gene)) + 
# print(ggplot(df, aes(width, color = gene)) + 
#   geom_density(alpha = 0.3) + 
#   scale_colour_Publication() +
#   theme_bw() + 
#   xlab("") + 
#   scale_x_continuous(limits = c(-50, 50)) + 
#   ggtitle("AFGHR" ) + 
#   theme(panel.grid.major = element_blank(),panel.grid.major.x = element_blank()) + 
#   theme(panel.grid.minor = element_blank()))
# dev.off()

#Plot the width accuracy
# df_list <- lapply(genes, function(g){
#   indel_events <- read.csv(file.path(input_dir,g,"indel_events.csv"),stringsAsFactors=FALSE)
#   if (nrow(indel_events) > 0){
#     qname_split <- strsplit(as.character(indel_events[,"read_id"]),":")
#     q <- do.call(rbind,qname_split)
#     indel_events <- cbind(q,indel_events[,4:9],stringsAsFactors=FALSE)
#     
#     mult_vec <- rep(1,nrow(indel_events))
#     mult_vec[indel_events[,9] == "D"] <- -1
#     indel_events[,8] <- as.numeric(indel_events[,8])
#     indel_events[,8] <- indel_events[,8] * mult_vec
#     
#     indel_events <- indel_events[abs(indel_events[,8]) == indel_events[,"width"],]
#     indel_events <- indel_events[indel_events[,9] == indel_events[,"indel_type"],]
#     accuracy <- table(factor(as.character(indel_events[,8]),levels=-max_width:max_width))
#     return(data.frame(indel=as.numeric(names(accuracy)),accuracy=as.numeric(accuracy),gene=g))
#   } else {
#     return(list())
#   }
# })
# df <- do.call(rbind,df_list)
# p <- ggplot(df, aes(x=indel, y=accuracy, group=gene,colour = gene))
# pdf(file.path(out_dir,"Width_Accuracy.pdf"),width=8,height=6)
# print(p + ggtitle("AGV28") + geom_line(size=1) + scale_colour_Publication() + theme_bw() + scale_x_continuous(breaks = seq(-150, 150, by = 20)))
# dev.off()

#Determine which reads are NA
# genes <- list.files(input_dir)
# df_list <- lapply(genes,function(g){
#   WT_info <- read.csv(file.path(input_dir,g,"all_reads_WT_info.csv"), stringsAsFactors=FALSE)
#   WT_info <- WT_info[order(WT_info[,4],decreasing=T),]
#   WT_info <- WT_info[!duplicated(WT_info[,1]),]
#   if (nrow(WT_info) > 0){
#     return(WT_info[is.na(WT_info[,"WT"]),])
#   } else {
#     return (list())
#   }
# })
# df <- do.call(rbind,df_list)


#Determine which reads are mistakes that are not NA
# genes <- list.files(input_dir)
# df_list <- lapply(genes,function(g){
#   WT_info <- read.csv(file.path(input_dir,g,"all_reads_WT_info.csv"), stringsAsFactors=FALSE)
#   if (nrow(WT_info) > 0){
#     qname_split <- strsplit(WT_info[,1],":")
#     q <- do.call(rbind,qname_split)
#     WT_info <- cbind(q,WT_info[,2:3],stringsAsFactors=FALSE)
#     
#     #Mistakes that are not NA
#     WT_info_mistake <- WT_info[!is.na(WT_info[,"WT"]),]
#     WT_info_mistake[,9] <- WT_info_mistake[,9] == "W"
#     WT_info_mistake <- WT_info_mistake[WT_info_mistake[,9] != WT_info_mistake[,"WT"],]
#     
#     return()
#   } else {
#     return (list())
#   }
# })
