library("Rsamtools")

args <- commandArgs(trailingOnly = TRUE)
#mutation_info_file <- args[1]
mutation_info_file <- "/Users/mburger/Data/crispr/mouse2/point_mutation_info.txt"
#gene_info_file <- args[2]
gene_info_file <- "/Users/mburger/Data/crispr/mouse2/mouse2_point_mut_genes.csv"
#bam_folder <- args[3]
bam_folder <- "/Users/mburger/Data/crispr/mouse2/Reads"
#min_q <- as.numeric(args[4])
min_q <- 20
#max_d <- as.numeric(args[5])
max_d <- 300000
#output_folder <- args[6]
output_folder <- "/Users/mburger/Data/crispr/mouse2"

output_folder <- file.path(output_folder,"point_mutation_quant")
mutect_folder <- file.path(output_folder,"mutect")
if (!file.exists(output_folder)){dir.create(output_folder)}
if (!file.exists(mutect_folder)){dir.create(mutect_folder)}

point_mut_genes <- read.csv("~/Data/crispr/mouse2/mouse2_point_mut_genes.csv", header=TRUE, row.names = 1, stringsAsFactors=FALSE)
mutation_info <- read.table(mutation_info_file,header=TRUE,stringsAsFactors=F,colClasses = "character",row.names=1)

all_samples <- unlist(strsplit(mutation_info[,"samples"],","))

output <- matrix(rep("",length(all_samples)*nrow(mutation_info)),length(all_samples),nrow(mutation_info))
colnames(output) <- rownames(mutation_info)
rownames(output) <- all_samples

for (g in rownames(mutation_info)){
  samples <- unlist(strsplit(mutation_info[g,"samples"],","))
  for (s in samples){
    refseq_interval <- paste0(point_mut_genes[g,"CHR"],":",point_mut_genes[g,"START"],"-",point_mut_genes[g,"END"])
    bam_file <- file.path(bam_folder,paste0(s,".sorted.bam"))
    cmd <- paste0("/System/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home/bin/java -Xmx2g -jar /Users/mburger/muTect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence /Users/mburger/Data/Homo_sapiens_assembly19.fasta --cosmic /Users/mburger/Data/hg19_cosmic_v54_120711.vcf --dbsnp /Users/mburger/Data/dbsnp_134_b37.leftAligned.vcf --intervals ",refseq_interval," --input_file:tumor ",bam_file," --enable_extended_output --out ", file.path(mutect_folder,paste0(s,".call_stats.txt"))," --vcf ",file.path(mutect_folder,paste0(s,".vcf"))," --coverage_file ",file.path(mutect_folder,paste0(s,".coverage.wig.txt"))) 
    system(cmd)
    
    position <- as.numeric(mutation_info[i,"start"])
    res <- pileup(file.path(bam_folder,paste0(s,".sorted.bam")),scanBamParam=ScanBamParam(which=GRanges(mutation_info[i,"chr"],IRanges(position,position))),pileupParam=PileupParam(max_depth=max_d,min_mapq=min_q,min_base_quality=min_q))
    total_reads <- sum(res$count)
    wt_reads <- res[res[,"nucleotide"] == mutation_info[i,"wt"],"count"]
    mut_reads <- res[res[,"nucleotide"] == mutation_info[i,"mut"],"count"]
    if (length(wt_reads) == 0){wt_reads <- 0}
    if (length(mut_reads) == 0){mut_reads <- 0}
    mut_frac <- round(mut_reads/total_reads,5)
    if (is.nan(mut_frac)){mut_frac <- 0}
    output[s,i] <- paste0(mut_frac,";",total_reads)
  }
}

write.csv(output,file=file.path(output_folder,paste0("pileup_q",min_q,".csv")),row.names=T)

