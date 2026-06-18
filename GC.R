library(data.table)
file <- fread("/mnt/SSD_4Tb/kmer/results/kmer_stats/joined_statistics0606_cleaned.tsv", header = T)
set.seed(123)
random_indices10 <- sample(1:nrow(file), 13422592)
random_sample10 <- file[random_indices10, ]
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

# Загрузка пакета Biostrings
library(Biostrings)
gc_content_biostrings <- function(sequence) {
  # Преобразование последовательности в объект DNAString
  dna_string <- DNAString(sequence)

  # Вычисление GC-состава
  gc_percentage <- letterFrequency(dna_string, "GC", as.prob = TRUE) * 100

  return(gc_percentage)
}

# Применяем функцию к столбцу k-меров и создаем новый столбец с GC-составом
random_sample10$gc_content <- sapply(random_sample10$kmer, gc_content_biostrings)
write.table(random_sample10, "../results/kmer_stats/subset10.txt", quote = FALSE)
