

split_matrix <- function(m,numSplits){
  mats <- list()
  for(j in 1:numSplits){
    mats[[j]] <- do.call(cbind, lapply(c(1:ncol(m)), function(i) as.numeric(matrix(unlist(strsplit(m[,i],";")), ncol=2, byrow=TRUE)[,j])))
    colnames(mats[[j]]) <- colnames(m)
    rownames(mats[[j]]) <- rownames(m)
  }  
  mats
}

combine_matrix <- function(m_list){
  m_1 <- m_list[[1]]
  for (i in 2:length(m_list)){
    for (j in 1:ncol(m_1)){
      m_1[,j] <- paste(m_1[,j],m_list[[i]][,j],sep=";")
    }
  }
  
  return(m_1)
}

individual_indel_frac <- function(results_dir,sample,gene){
  
  indel_mat <- get_read_level_wt_call(results_dir,sample,gene,"WT_info")
  
  if (is.null(indel_mat)){
    return(NA)
  } else if (nrow(indel_mat) == 0){
    return(NA)
  }
  
  indel_mat <- indel_mat[order(indel_mat[,"WT"]),]
  indel_mat <- indel_mat[!duplicated(indel_mat[,"read_id"]),]
  calls <- indel_mat[!is.na(indel_mat[,"WT"]),"WT"]
  mut <- sum(calls == 0)
  wt <- sum(calls == 1)
  total <- mut + wt
  if (mut > 0){
    frac <- mut / total
  } else {
    frac <- 0
  }
  frac <- paste0(signif(frac, digits = 4),";",total)
  return(frac)
}


genes_per_sample_indel_frac <- function(results_dir,sample,gene_list){
  
  sample_row <- sapply(gene_list, function(g) individual_indel_frac(results_dir,sample,g)) 
  return(sample_row)
  
}


create_indel_frac_table <- function(results_dir,sample_list,gene_list){
  
  indel_frac <- sapply(sample_list, function(s) genes_per_sample_indel_frac(results_dir,s,gene_list))
  return(t(indel_frac))
}

negative_exact_test <- function(indel_fractions,indel_totals,neg_samples){
  pos_exp <- round(indel_fractions * indel_totals, digits = 0)
  neg_exp <- round((1 - indel_fractions) * indel_totals, digits = 0)
  
  #pos_cntrl <- ceiling(colMeans(pos_exp[neg_samples,]))
  #neg_cntrl <- ceiling(colMeans(neg_exp[neg_samples,]))
  pos_cntrl <- rep(NA,ncol(indel_fractions))
  neg_cntrl <- rep(NA,ncol(indel_fractions))
  for (i in 1:ncol(pos_exp)){
    gene <- colnames(indel_fractions)[i]
    neg_samples_gene <- rownames(neg_samples)[neg_samples[,gene] == 1]
    pos_cntrl[i] <- ceiling(mean(pos_exp[neg_samples_gene,gene]))
    neg_cntrl[i] <- ceiling(mean(neg_exp[neg_samples_gene,gene]))
  }
  
  pvals <- pos_exp
  
  for (i in 1:nrow(pos_exp)){
    for (j in 1:ncol(pos_exp)){
      fisher_matrix <- cbind(c(pos_cntrl[j],neg_cntrl[j]),c(pos_exp[i,j],neg_exp[i,j]))
      if(any(is.nan(fisher_matrix))){
        pvals[i,j] <- NA
      } else {
        pvals[i,j] <- fisher.test(fisher_matrix, alternative="less")$p.value
      }
    }
  }
  
  qvals <- p.adjust(unlist(pvals),method="BH")
  qvals <- matrix(qvals,nrow=dim(pvals)[1],ncol=dim(pvals)[2],byrow=F)
  rownames(qvals) <- rownames(pvals)
  colnames(qvals) <- colnames(pvals)
  
  return(qvals)
}