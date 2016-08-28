
args <- commandArgs(trailingOnly = TRUE)
code_dir <- args[1]
#code_dir <- "/Users/mburger/Documents/DNAseq_pipelines/targeted_CRISPR_SE"
#setwd(code_dir)
source(file.path(code_dir,"get_data.R"))
source(file.path(code_dir,"LOF.R"))
source(file.path(code_dir,"load_annotations.R"))
source(file.path(code_dir,"indel_frac.R"))
source(file.path(code_dir,"indel_descriptions.R"))


project_dir <- args[2]
#id <- "AL23V"
#base_dir <- "/Users/mburger/Data/crispr"
#project_dir <- file.path(base_dir,id)
pysam_outdir <- file.path(project_dir,"VariantCalls")
cutsite_annot <- args[3]
#cutsite_annot <- file.path(base_dir,id,paste0(id,"_indel_cut_interval.csv"))
neg_cntrl <- args[4]
#neg_cntrl <- file.path(base_dir,id,paste0(id,"_negative_controls.csv"))

# Import intervals from csv
target_dict <- load_intervals(cutsite_annot)

# Get sample and gene names
sample_names <- list.files(pysam_outdir)
gene_names <- names(target_dict)

# Create output dirs
base_outdir <- file.path(project_dir,"indel_quant")
dir.create(base_outdir, recursive=T)

descriptions_outdir <- file.path(base_outdir,"descriptions")
tables_outdir <- file.path(base_outdir,"tables")
plots_outdir <- file.path(base_outdir, "plots_trial")
dir.create(plots_outdir, recursive=T)
dir.create(descriptions_outdir, recursive=T)
dir.create(tables_outdir, recursive=T)

indel_frac_table <- create_indel_frac_table(pysam_outdir,sample_names,gene_names)
indel_frac_table[is.na(indel_frac_table)] <- "0;0"
write.csv(indel_frac_table,file=file.path(tables_outdir,"raw_indel_fractions_totals.csv"),row.names=T)

indel_frac_split <- split_matrix(indel_frac_table,2)
indel_fractions <- indel_frac_split[[1]]
indel_totals <- indel_frac_split[[2]]
write.csv(indel_totals,file=file.path(tables_outdir,"raw_indel_totals.csv"),row.names=T)

if (neg_cntrl != "NA"){
  # Get negative control sample names
  #neg_samples <- read.csv(neg_cntrl, quote="\"", comment.char="", stringsAsFactors=FALSE)[,1]
  neg_samples <- read.csv(neg_cntrl, row.names=1, stringsAsFactors=FALSE)
  rownames(neg_samples) <- gsub(" ","",rownames(neg_samples))
  neg_samples <- neg_samples[rownames(neg_samples) %in% rownames(indel_frac_table),]
#   if (typeof(neg_samples) != "character"){
#     neg_samples <- as.character(neg_samples)
#   }
  
  qvals <- negative_exact_test(indel_fractions,indel_totals,neg_samples)
  write.csv(qvals,file=file.path(tables_outdir,"q_values.csv"),row.names=T)
  
  indel_fractions[qvals > .05] <- 0
  indel_fractions[is.na(qvals)] <- NA
  write.csv(indel_fractions,file=file.path(tables_outdir,"sig_indel_fractions.csv"),row.names=T)
  
  indel_frac_split[[1]] <- indel_fractions
  sig_indel_frac_totals <- combine_matrix(indel_frac_split)
  write.csv(sig_indel_frac_totals,file=file.path(tables_outdir,"sig_indel_fractions_totals.csv"),row.names=T)
  
} else {
  write.csv(indel_fractions,file=file.path(tables_outdir,"raw_indel_fractions.csv"),row.names=T)
}



create_indel_descriptions(pysam_outdir,descriptions_outdir,sample_names,gene_names)
