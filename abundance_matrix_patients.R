library(data.table)
library(tidyverse)
library(purrr)
library(tibble)
library(dplyr)
library("stringr")

start.time <- Sys.time()

#start_row<-1
#end_row<-10000000
#total_rows<-end_row-start_row+1
header <- fread("../files/result_abundance_1_10.txt", nrows = 1, fill = TRUE)
#print(header)


data <- fread("../files/result_abundance_1_10.txt",skip = 1)
setnames(data, colnames(header))
meta_data_new_kyh <-read.csv("../wgs_julia.csv")
grouped_data <- meta_data_new_kyh %>%
  group_by(ID_KYH) %>%
  summarise(Values = paste(sample_id, collapse = ";")) %>%
  ungroup()


select_and_rename_columns <- function(id, columns, df) {
  columns_to_select <- str_split(columns, ";")[[1]]
  selected_data <- df %>% select(all_of(columns_to_select))
  selected_data[is.na(selected_data)] <- 0
  colnames(selected_data) <- paste0(id, "_", colnames(selected_data))
  return(selected_data)
}



result_list <- map2(grouped_data$ID_KYH, grouped_data$Values, ~select_and_rename_columns(.x, .y, df = data))
library(purrr)
library(readr)


update_column_names <- function(df, name_file, sign) {

  required_names <- read_lines(name_file)
  

  new_colnames <- sapply(colnames(df), function(colname) {

    parts <- str_split(colname, "_")[[1]]
    parts2 <- paste0(parts[2], "_", parts[3])
   
    if (length(parts) > 1 && parts2 %in% required_names) {
      return(paste0(colname, sign))
    } else {
      return(colname)
    }
  })

  colnames(df) <- new_colnames
  return(df)
}

calculate_fold_change <- function(df) {
  colnames_list <- colnames(df)[-1] 
  

  parts <- str_split(colnames_list, "_")
  first_parts <- sapply(strsplit(colnames_list, "_"), `[`, 1)
  second_parts <- paste0(sapply(strsplit(colnames_list, "_"), `[`, 2), "_", sapply(strsplit(colnames_list, "_"), `[`, 3))
  last_parts <- sapply(strsplit(colnames_list, "_"), `[`, 4)
  parts_df <- data.frame(do.call(rbind, parts), stringsAsFactors = FALSE)
  colnames(parts_df) <- c("Part1", "Part2", "Part3", "Year")
  

  filtered_df <- parts_df %>% filter(Year %in% c("2017", "2022"))
  

  unique_parts <- unique(filtered_df$Part1)
  
  
  result <- df %>% select(Rownames)
  
  for (part in unique_parts) {
  
    cols_2017 <- colnames_list[filtered_df$Part1 == part & parts_df$Year == "2017"]
    cols_2022 <- colnames_list[filtered_df$Part1 == part & parts_df$Year == "2022"]
    if (length(cols_2017) > 0 && length(cols_2022) > 0) {
      for (col_2017 in cols_2017) {
        col_2022 <- gsub("2017", "2022", col_2017)
        col_2022_first <- sapply(strsplit(col_2022, "_"), `[`, 1)
        col_2022_last <- sapply(strsplit(col_2022, "_"), `[`, 4)
        
        if (col_2022_first %in% first_parts) {
          index2 <- which(first_parts == col_2022_first)
         
          if (length(index2) > 0) {
            if (last_parts[index2[1]] == "2022") {
              index <- index2[1]
            } else {
              if (length(index2) > 1) {
                index <- index2[2]
              }
            }
          }
          col_2022_1 <- paste0(first_parts[index], "_", second_parts[index], "_", last_parts[index])
          
          num1 <- as.numeric(df[[col_2022_1]])
          num2 <- as.numeric(df[[col_2017]])
          
          fold_change_col <- num1 / (num2 + 1)
          mean_col <- (num1 + num2) / 2
          
          if (length(fold_change_col) > 0 && !any(is.na(fold_change_col))) {
            a <- paste(part, "fold_change", sep = "_")
            result[[a]] <- fold_change_col
          }
          if (length(mean_col) > 0 && !any(is.na(mean_col))) {
            b <- paste(part, "mean_col", sep = "_")
            result[[b]] <- mean_col
          }
        }
      }
    }
  }
  
  return(result)
}

process_element <- function(i, result_list, data) {
  result_df <- data.frame(result_list[[i]])
  result_df$Rownames <- data$Rownames
  
  df_combined <- result_df %>%
    group_by(Rownames) %>%
    summarise(across(everything(), ~paste(na.omit(.), collapse = ","), .names = "{col}")) %>%
    ungroup()
  names(df_combined) <- sub("^X", "", names(df_combined))
  print(df_combined)
  
  df_combined <- update_column_names(df_combined, "../files/names_2017.txt", "_2017")
  df_combined <- update_column_names(df_combined, "../files/names_2022.txt", "_2022")
  
  df <- calculate_fold_change(df_combined)
  
  return(df)
}
result_dfs <- vector("list", length(result_list))



result_dfs <- lapply(1:length(result_list), process_element, result_list = result_list, data = data)
df <- Reduce(function(x, y) merge(x, y, by = "Rownames", all = TRUE), result_dfs)
write.table(df,"../files/result_abundance_1_10_patients.txt",quote = FALSE)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
print(time.taken)
