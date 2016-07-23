
##
# Compiles the insertion and deletion lists for a given sample and gene.
# Returns null if there are not any files matching this description 
#
get_sample_gene_vcs <- function(results_dir,sample,gene){
  
  #check dir format
  #results_dir <- ifelse(substr(results_dir, nchar(results_dir),nchar(results_dir)) == "/", paste(results_dir,sample,"/",gene,"/",sep=""), paste(results_dir,"/",sample,"/",gene,"/",sep=""))
  
  file_list <- list.files(results_dir,pattern="_indels.csv")
 
  vcfs <- lapply(file_list,function(f) read.csv(file.path(results_dir,f),stringsAsFactors=F))
  vcfs <- do.call(rbind,vcfs)
}

get_read_level_indels <- function(input_dir,sample,gene,fname){
  
  #check dir format and indel_type
  #file_dir <- ifelse (substr(input_dir, nchar(input_dir),nchar(input_dir)) == "/",paste(input_dir,sample,"/",gene,"/",sep=""),paste(input_dir,"/",sample,"/",gene,"/",sep=""))
  file_dir <- file.path(input_dir,sample,gene)
  file_list <- list.files(file_dir,pattern=paste0(fname,".csv"))
  
  vcfs <- lapply(file_list,function(f) read.csv(file.path(file_dir,f)))
  vcfs <- do.call(rbind,vcfs)
  return(vcfs)
}

get_read_level_wt_call <- function(input_dir,sample,gene,fname){
  
  #check dir format and indel_type
  #file_dir <- ifelse (substr(input_dir, nchar(input_dir),nchar(input_dir)) == "/",paste(input_dir,sample,"/",gene,"/",sep=""),paste(input_dir,"/",sample,"/",gene,"/",sep=""))
  file_dir <- file.path(input_dir,sample,gene)
  
  file_list <- list.files(file_dir,pattern=paste0(fname,".csv"))
  
  vcfs <- lapply(file_list,function(f) read.csv(file.path(file_dir,f)))
  vcfs <- do.call(rbind,vcfs)
  return(vcfs)
}