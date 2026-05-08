# ============================================================================
# Natural Connectivity Calculation for Network Robustness Analysis
# Author: Yujie Wang
# Date: 2026-05-08
# 
# Description:
#   This script calculates natural connectivity (a spectral measure of network 
#   robustness) for a microbial co-occurrence network. It computes the natural 
#   connectivity of the original network and after sequential random node removal 
#   at specified proportions (10%, 20%, 30%, 40%, 50%). Natural connectivity 
#   is defined as the average of the eigenvalues of the adjacency matrix.
#
# Parameters:
#   input_file: "227edge_RH_3h.csv" - Network edge list in CSV format
#                Columns: Source, Target, Weight (and other attributes)
#   output_file: "natural_connectivity_results_RH_3h.csv" - Results output file
#   removal_proportions: c(0.1, 0.2, 0.3, 0.4, 0.5) - Fractions of nodes to remove
#   random_seed: Not explicitly set (uses default R random number generator)
#
# Output:
#   Data frame with two columns:
#     - ProportionRemoved: Fraction of nodes removed (0.1 to 0.5)
#     - NaturalConnectivity: Natural connectivity value after node removal
#
# Required Packages:
#   igraph
# 
# Reference:
#   Natural connectivity as a measure of network robustness.
#   Jun, W. et al. (2010). Physica A: Statistical Mechanics and its Applications.
# ============================================================================


# Load required library
library(igraph)

# Read CSV file containing network edges
# Modify the file path as needed for your environment
file_path <- "227edge_RH_3h.csv"  
df <- read.csv(file_path, stringsAsFactors = FALSE)

# Create undirected graph from edge list
G <- graph_from_data_frame(d = df, directed = FALSE)

# Function to calculate natural connectivity
# Natural connectivity = average of eigenvalues of the adjacency matrix
natural_connectivity <- function(G) {
  # Get adjacency matrix of the graph
  A <- as_adjacency_matrix(G, sparse = FALSE)
  
  # Calculate eigenvalues of the adjacency matrix
  eigenvalues <- eigen(A)$values
  
  # Ensure eigenvalues are non-negative
  eigenvalues <- Re(eigenvalues)  # Extract real part only
  eigenvalues <- eigenvalues[eigenvalues > 0]  # Remove negative and zero values
  
  # Calculate natural connectivity
  n <- vcount(G)  # Number of nodes in the graph
  if (n == 0) {
    return(0.0)
  }
  phi <- sum(eigenvalues) / n
  return(phi)
}

# Define node removal proportions for robustness testing
# Test robustness at 10%, 20%, 30%, 40%, and 50% node removal
proportions <- c(0.1, 0.2, 0.3, 0.4, 0.5)  # 10%, 20%, 30%, 40%, 50%
results <- data.frame(ProportionRemoved = proportions, NaturalConnectivity = numeric(length(proportions)))

# Calculate natural connectivity after removing specified proportion of nodes
for (i in seq_along(proportions)) {
  p <- proportions[i]
  # Randomly remove specified proportion of nodes
  nodes_to_remove <- sample(V(G), size = round(p * vcount(G)))
  G_reduced <- delete_vertices(G, nodes_to_remove)
  results$NaturalConnectivity[i] <- natural_connectivity(G_reduced)
}

# Display results in console
print(results)

# Save results to CSV file
output_file_path <- "natural_connectivity_results_RH_3h.csv"  
write.csv(results, file = output_file_path, row.names = FALSE)
