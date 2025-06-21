# Load necessary libraries
library(igraph)
library(tidyverse)

# Load the data
# IMPORTANT! Screen your dataset first so that only relations with 
# P < 0.05 and r > 0.6 are retained for analysis.
data <- read_csv("<path-to-data.csv>")

# Create the graph
g <- graph_from_data_frame(data %>% select(<source-column>, <target-column>), directed = FALSE)

# Calculate metrics
betweenness_values <- betweenness(g)
closeness_values <- closeness(g)
degree_values <- degree(g)
mean_degree <- mean(degree_values)

# Create a data frame with calculated metrics
otu_metrics <- data.frame(OTU = V(g)$name,
                          Betweenness = betweenness_values,
                          Closeness = closeness_values,
                          Degree = degree_values,
                          Mean_Degree = mean_degree)

# Define top 30% and bottom 30% thresholds
top_30_mean_degree <- quantile(otu_metrics$Mean_Degree, 0.70)
top_30_degree_centrality <- quantile(otu_metrics$Degree, 0.70)
top_30_betweenness <- quantile(otu_metrics$Betweenness, 0.70)
top_30_closeness <- quantile(otu_metrics$Closeness, 0.70)
bottom_30_betweenness <- quantile(otu_metrics$Betweenness, 0.30)

# Identify hub OTUs (top 30% in degree centrality and betweenness)
hub_OTUs <- otu_metrics %>%
  filter(Degree >= top_30_degree_centrality & Betweenness >= top_30_betweenness)

# Identify keystone OTUs (top 30% in mean degree, top 30% in closeness, and bottom 30% in betweenness)
keystone_OTUs <- otu_metrics %>%
  filter(Mean_Degree >= top_30_mean_degree & 
           Closeness >= top_30_closeness & 
           Betweenness <= bottom_30_betweenness)