library(data.table)
library(tidyverse)
library(purrr)
library(tibble)
library(dplyr)

# Load necessary library
library(stringr)

start.time <- Sys.time()
header <- fread("../files/matrix_all_1_10.txt", nrows = 1)
#header <- fread("../matrix_100_test.txt",nrows = 1)
data <- fread("../files/matrix_all_1_10.txt",skip = 1)
#data <- fread("../matrix_100_test.txt",skip = 1)
# Присваиваем заголовок считанным строкам
setnames(data, colnames(header))

# Read the file
file_path <- "count_reads.txt"  # Change this to your actual file path
lines <- readLines(file_path)

# Extract basename and reads
basenames <- str_extract(lines, "(?<=/)[^/]+(?=.fastq.gz)")
reads <- as.integer(str_extract(lines, "\\d+(?= reads)"))

# Create dataframe
df <- data.frame(basename = basenames, reads = reads)
df2_filtered <- df[str_detect(df$basename, "_R1$"), ]
# Print dataframe
merge_cols <- colnames(data)[2:length(colnames(data))]  # Exclude the first column

transposed_data_df <- as.data.frame(t(data))
transposed_data_df <-transposed_data_df[2:nrow(transposed_data_df),]
names(transposed_data_df) <- data$RowNames
transposed_data_df$basename <-rownames(transposed_data_df)
df2_filtered$basename <- gsub("_R1$", "", df2_filtered$basename)

#transposed_data_df <- transposed_data_df[order(transposed_data_df$basename), ]
merged_df<-merge(x=transposed_data_df, y = df2_filtered, by = "basename", all = TRUE) #порядок в колонке дб тот же

last_column_values <- as.numeric(merged_df[, ncol(merged_df)])


matrix_values <- as.matrix(sapply(merged_df[, 2:(ncol(merged_df) - 1)], as.numeric))
last_values <- as.numeric(merged_df[, ncol(merged_df)])

# Perform the division using vectorized operations
result_matrix <- sweep(matrix_values, 1, last_values, "/")

merged_dt <- as.data.table(merged_df)
results_dt <- as.data.table(result_matrix)

# Bind the first column of merged_df with results
results <- cbind(merged_dt[, 1, with = FALSE], results_dt)
results_t <- t(results)
results_t <-as.data.frame(results_t)
colnames(results_t) <- as.character(results_t[1, ])

# Remove the first row
results_t <- results_t[-1, ]
write.table(results_t,"../files/result_abundance_1_10.txt", quote = FALSE)
#write.table(results_t,"../result_abundance_test.txt", quote = FALSE)
