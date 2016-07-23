##
# Takes a 4 column csv: [GENE,CHR,START,END]
# Returns a list indexed by gene name with value equal to a num vector of length 2, 
# which are start and end positions respectively 
#
load_intervals <- function(input_csv){
  
  if (!file.exists(input_csv)){print("Error: Target intervals csv annotation is missing!")}
  
  raw_table <- read.csv(input_csv,header=F,stringsAsFactors=F)
  colnames(raw_table) <- c("gene","chr","start","end")
  rownames(raw_table) <- raw_table[,1]
  target_dict <- lapply(1:nrow(raw_table),function(g) raw_table[g,3:4])
  names(target_dict) <- raw_table[,1]
  
  return(target_dict)
}