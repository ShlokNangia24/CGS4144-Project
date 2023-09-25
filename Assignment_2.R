data_dir <- file.path("R/CGS_data/data","SRP073813")
data_file <- file.path(data_dir,"SRP073813_Entrez_IDs.tsv")
metadata_file <- file.path(data_dir,"metadata_SRP073813.tsv")
df <- readr::read_tsv(data_file)
df_meta <- readr::read_tsv(metadata_file)
# Assuming df is your loaded data and removing non-numeric columns
# Load the required libraries
library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::select)
drops <- c("Ensembl","all_entrez_ids")
df <- df[ ,!(names(df) %in% drops)]

# Remove non-numeric columns from df
df_numeric <- df %>% select(where(is.numeric))

# Log-scale the data (adding a small constant to avoid log(0))
log_df <- log2(df_numeric + 1)

# Calculate per-gene median expression ranges
ranges <- apply(log_df, 1, function(row) {
  median_range <- diff(range(row))
  return(median_range)
})


# Make a density plot
plot <- ggplot(data.frame(ranges), aes(x=ranges)) +
  geom_density(fill="blue", alpha=0.5) +
  labs(title="Density Plot of Per-Gene Median Expression Ranges",
       x="Median Expression Range", y="Density")

print(plot)

# Part two ----------------------------------------------------------------
cts <- as.matrix(df)
coldata <- readr::read_tsv(metadata_file)
coldata <- coldata[,c("refinebio_subject","refinebio_accession_code")]
coldata$refinebio_subject <- factor(coldata$refinebio_subject)
coldata$refinebio_accession_code <- factor(coldata$refinebio_accession_code)
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
vsd <- vst(dds, blind=FALSE)
library("vsn")
plotPCA(vsd, intgroup=c("condition", "type"))


