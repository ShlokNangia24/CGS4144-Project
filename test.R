# Load necessary libraries
library(ggalluvial)
library(ggplot2)
data_dir <- file.path("R/CGS_data/data","SRP073813")
data_file <- file.path(data_dir,"Sankey_plot.csv") 
# Read the data
data <- read.csv(data_file, row.names = 1)

# Melt the data for plotting
data_melted <- reshape2::melt(data)
numGenes <- rownames(data_melted)
# Plot the alluvial diagram
ggplot(data_melted, aes(axis1 = variable, axis2 = as.factor(value), y = numGenes)) +
  geom_flow(aes(fill = as.factor(value))) +  # Use distinct colors for the flows
  geom_stratum(aes(fill = variable)) +       # Add colors to the strata
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  # Add labels to the strata
  theme_minimal() +
  labs(title = "Alluvial diagram of cluster memberships across different gene subsets",
       fill = "Cluster") +  # Add legend title
  theme(legend.position = "bottom")  # Position the legend at the bottom
