args <- commandArgs(trailingOnly = TRUE)

output_folder <- args[1]
config_file <- args[2]
picard_path <- args[3]

#known_indel_quant <- file.path(output_folder,"known_indel_quant")
#ref_folder <- file.path(known_indel_quant,"ref")
ref_folder <- output_folder
dir.create(ref_folder,recursive=T)
config_name <- gsub(".txt","",basename(config_file))
out_file <- file.path(ref_folder,paste0(config_name,"_customRef.fasta"))

sequences <- read.delim(config_file, header=FALSE, stringsAsFactors=FALSE)

lengths <- nchar(sequences[,4])
sequences <- cbind(sequences,lengths)

rownames(sequences) <- paste0(">",sequences[,3]," dna:chromosome chromosome:GRCh37:",sequences[,1],":",sequences[,2],":",sequences[,2] + sequences[,5],":1")

for (i in 1:nrow(sequences)){
  write(rownames(sequences)[i],file=out_file,append=T)
  write(sequences[i,4],file=out_file,append=T)
}

cmd <- paste0("java -jar ",picard_path," CreateSequenceDictionary REFERENCE=",out_file," OUTPUT=",file.path(ref_folder,"custom_ref.dict")," GENOME_ASSEMBLY=GRCh37")
system(cmd)

cmd <- paste0("samtools faidx ",out_file)
system(cmd)

cmd <- paste0("bwa index ",out_file)
system(cmd)
