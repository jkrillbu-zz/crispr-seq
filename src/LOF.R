

##
# Calculates a single LOF fraction treating each read individually, not going per unique indel
# Depends on the LOF column being a logical and the count column being a num
# Returns zero if there aren't any indel files
#
individual_LOF_frac <- function(results_dir,sample,gene){
  
  indel_mat <- get_sample_gene_vcs(results_dir,sample,gene)
  
  if (is.null(indel_mat)){
    return(0)
  }
  
  #frac <- sum(indel_mat[indel_mat[,"LOF"],"count"]) / sum(indel_mat[,"count"])
  frac <- sum(indel_mat[,"LOF"]) / nrow(indel_mat)
  #frac <- nrow(indel_mat)
  
  return(frac)
}

##
# Calls the individual_LOF_frac function for each gene per sample
# Returns a gene named num vector
#
genes_per_sample_LOF_frac <- function(results_dir,sample,gene_list){
  
  sample_row <- sapply(gene_list, function(g) individual_LOF_frac(results_dir,sample,g)) 
  return(sample_row)
  
}

##
# Calls the sample row generating LOF frac function for each sample
#
create_LOF_frac_table <- function(results_dir,sample_list,gene_list){
  
  LOF_frac <- sapply(sample_list, function(s) genes_per_sample_LOF_frac(results_dir,s,gene_list))
  return(t(LOF_frac))
}