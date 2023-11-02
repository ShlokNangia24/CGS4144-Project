
# Part 1 ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)

data_dir <- file.path("R/CGS_data/data","SRP073813")
data_file <- file.path(data_dir,"nacc_Var.csv") 
metadata_file <- file.path(data_dir,"nacc_meta.csv")
meta_data <- read.csv(metadata_file)
df <- read.csv(data_file, row.names = "X")
df <- na.omit(df)
df <- head(df,5000)
dist_matrix <- dist(df, method="euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, main="Hierarchical Clustering Dendrogram", labels = FALSE )
heatmap(as.matrix(df))

annotation_df <- meta_data %>%
  # select only the columns we need for annotation
  dplyr::select(
    ID,
    refinebio_subject
  )
groupCluster <- annotation_df
groupCluster$refinebio_subject <-replace(groupCluster$refinebio_subject,  groupCluster$Group == "nacc_control", 1)
groupCluster$refinebio_subject <-replace(groupCluster$refinebio_subject,  groupCluster$Group == "nacc_schizophrenia", 2)
groupCluster$refinebio_subject <-replace(groupCluster$refinebio_subject,  groupCluster$Group == "nacc_major depression", 3)
groupCluster$refinebio_subject <-replace(groupCluster$refinebio_subject,  groupCluster$Group == "nacc_bipolar disorder", 4)
originalGroups <- as.matrix(groupCluster[1:114, 2])
rownames(originalGroups) <- groupCluster$ID

temp<- cutree(hc, k =3)
chisq1 <- chisq.test(df)
p = c(chisq1$p.value)
pvalues = c(chisq1$p.value)
adjusted <-p.adjust(pvalues, method = "fdr", n = length(p))
chisq1a <- chisq1
chisq1a$p.value <- adjusted[2]

adj_pValues <- data.frame(
  regular = c(chisq1$p.value),
  adjusted = c(chisq1a$p.value)
)
adj_pValues
chisq1
