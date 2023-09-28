# Part One ----------------------------------------------------------------

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
df <- na.omit(df)

# Remove non-numeric columns from df
df_numeric <- df %>% select(where(is.numeric))
# Log-scale the data (adding a small constant to avoid log(0))
df_numeric <- (df_numeric - min(df_numeric) +2)
log_df <- log2(df_numeric)

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

# Part five ----------------------------------------------------------------
data_file <- file.path(data_dir,"nacc_log.csv")
metadata_file <- file.path(data_dir,"nacc_meta.csv")
df <- read.csv(data_file)
df_meta <- read.csv(metadata_file)
library(gprofiler2)
gostres <- gost(query = df$X, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
gostplot(gostres, capped = FALSE, interactive =FALSE)

# Part Six ----------------------------------------------------------------
BiocManager::install("GOplot")
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(gene          = gene,
                universe      = df$X,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego,3)
goplot(ego)
