start.time <- Sys.time()
library(data.table)
library(dplyr)
library(tools)

# List of file names (adjust the pattern if needed)
directory <- "/mnt/16TB/PROJECTS/arhangels_kmer/snake/results_16mer"
# Получение списка всех файлов CSV в директории
files <- list.files(directory, pattern = "\\.kmer$", full.names = TRUE)
#print(paste0("files ",files)

merged_df <- fread(files[1])
count_column_index0 <- 2
kmer_column_index0<- 1  
file_name0 <-file_path_sans_ext(basename(files[1]))
 

colnames(merged_df)[count_column_index0] <- file_name0
colnames(merged_df)[kmer_column_index0] <- "kmer"
 #dfs <- lapply(files, fread)
#print(paste0("dfs ",dfs))
# Use file names without extensions as column names for the second column of each data frame

for (file in files[-1]) {
  temp_df <- fread(file)
  file_name <-file_path_sans_ext(basename(file))
  count_column_index <- 2
  kmer_column_index<- 1
  colnames(temp_df)[count_column_index] <- file_name
  colnames(temp_df)[kmer_column_index] <- "kmer"
  merged_df <- full_join(merged_df, temp_df, by = names(merged_df)[1])
}


#for (i in seq_along(dfs)) {
  # Extract file name without extension
#  file_name <- file_path_sans_ext(basename(files[i]))
 # print(file_name)
  # Rename the second column to the file name
#  count_column_index <- 2
#  kmer_column_index<- 1
#  colnames(dfs[[i]])[count_column_index] <- file_name
 # colnames(dfs[[i]])[kmer_column_index] <- "kmer"
#
#}

# Initialize merged_df with the first data frame
#merged_df <- dfs[[1]]

# Perform a full join for each subsequent data frame
#for (i in 2:length(dfs)) {
 # merged_df <- full_join(merged_df, dfs[[i]], by = "kmer")
#}

# Print the final merged data frame
# Сохранение результата в новый CSV файл
fwrite(merged_df, "../results/combined_R_16.csv",sep = "\t")
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

