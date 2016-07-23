library("ShortRead")

##
# Takes a matrix and filters it to contain only rows that overlap the given interval
# @param interval_start : int 
# @param interval_end: int 
# @param data_mat: matrix or dataframe to filter
# @param data_start_col: column name or number of start positions in the data_mat
# @param data_end_col: column name or number of end positions in the data_mat
#
filter_by_interval <- function(interval_start,interval_end,data_mat,data_start_col,data_end_col){
  
  event_ranges <- IRanges(as.numeric(data_mat[,data_start_col]),as.numeric(data_mat[,data_end_col]))
  event_tree <- IntervalTree(event_ranges)
  
  query <- IRanges(as.numeric(interval_start),as.numeric(interval_end))
  f <- findOverlaps(query,event_tree,type="any",select="all")
  
  return(data_mat[f@subjectHits,])
  
}

individual_description <- function(input_dir,output_dir,sample,gene){
  
  combined_raw_calls <- get_read_level_indels(input_dir,sample,gene,"indel_events")
  #filtered_raw_calls <- filter_by_interval(target_intervals[[gene]]["start"],target_intervals[[gene]]["end"],combined_raw_calls,"start","end")
  filtered_raw_calls <- combined_raw_calls[combined_raw_calls[,"onTarget"] == 1,]
  
  final_header <- c("id","indel_type","start","end","width","count","LOF","indel_desc")
  
  if (nrow(filtered_raw_calls) > 0){
  
    #Give each indel an ID 
    id <- paste0(filtered_raw_calls[,"indel_type"],":",filtered_raw_calls[,"start"],":",filtered_raw_calls[,"end"])
    
    #Count occurances of IDs
    count <- table(id)
    count <- count[order(count,decreasing=T)]
    
    #Format columns for output
    spl <- strsplit(names(count),":")
    names(spl) <- names(count)
    df <- t(data.frame(spl))
    colnames(df) <- c("indel_type","start","end")
    width <- as.numeric(df[,"end"]) - as.numeric(df[,"start"]) + 1
    LOF <- sapply(width, function(x){if (x > 15 || (x %% 3) != 0){return(TRUE)} else {return(FALSE)}}) 
    df <- cbind(df,width,count,LOF) 
    indel_desc <- paste(df[,"indel_type"],":",df[,"width"],sep="")
    id <- paste(indel_desc," (",df[,"start"],")",sep="")
    df <- cbind(id,df,indel_desc)
  } else {
    df <- read.table(textConnection(""), col.names = final_header)
  }
  
  if (colnames(df) != final_header){
    stop("Indel description header output does not match")
  }
  dir.create(output_dir, recursive=T)
  write.table(df,file=file.path(output_dir,paste0(sample,"_",gene,"_indels.csv")),sep=",",row.names=F,col.names=T)

}

##
# Calls description function on each gene
#
per_sample_descriptions <- function(input_dir,output_dir,sample,gene_list){
  
  dir.create(output_dir, recursive=T)
  invisible(lapply(gene_list, function(g) individual_description(input_dir,file.path(output_dir,g),sample,g)))
  
}

##
# Calls description function on each sample
#
create_indel_descriptions <- function(input_dir,output_dir,sample_list,gene_list){
  
  invisible(lapply(sample_list, function(s) per_sample_descriptions(input_dir,file.path(output_dir,s),s,gene_list)))
  
}

