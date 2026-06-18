start.time <- Sys.time()
library(data.table)
library(dplyr)
library(tools)


directory <- "/mnt/16TB/PROJECTS/arhangels_kmer/snake/results_16mer"

files <- list.files(directory, pattern = "\\.kmer$", full.names = TRUE)


merged_df <- fread(files[1])
count_column_index0 <- 2
kmer_column_index0<- 1  
file_name0 <-file_path_sans_ext(basename(files[1]))
 

colnames(merged_df)[count_column_index0] <- file_name0
colnames(merged_df)[kmer_column_index0] <- "kmer"


for (file in files[-1]) {
  temp_df <- fread(file)
  file_name <-file_path_sans_ext(basename(file))
  count_column_index <- 2
  kmer_column_index<- 1
  colnames(temp_df)[count_column_index] <- file_name
  colnames(temp_df)[kmer_column_index] <- "kmer"
  merged_df <- full_join(merged_df, temp_df, by = names(merged_df)[1])
}



fwrite(merged_df, "../results/combined_R_16.csv",sep = "\t")
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

