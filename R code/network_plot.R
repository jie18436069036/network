# ============================================================================
# Network Construction and Edge/Node File Generation
# Author: Yujie Wang
# Date: 2026-05-08
# 
# Description:
#   This script constructs a microbial co-occurrence network based on Spearman 
#   correlation analysis. It filters OTUs/ASVs at the genus level, calculates 
#   pairwise correlations, applies significance thresholds, and generates 
#   edge and node files compatible with Gephi and other network visualization 
#   tools.
#
# Parameters:
#   input_otu_file: "network_RH_3h_1.txt" - OTU/ASV abundance table 
#                   (rows = taxa, columns = samples)
#   input_taxonomy_file: "Taxonomy.txt" - Taxonomy/attribute file for nodes
#   output_edge_file: "221edge_RH_3h.csv" - Network edge list with weights
#   output_node_file: "221node_with_attributes_RH_3h.csv" - Node list with attributes
#   correlation_method: "spearman" - Correlation method
#   p_value_threshold: 0.05 - Significance threshold for p-values
#   correlation_threshold: 0.75 - Absolute correlation coefficient cutoff
#   p_adjust_method: "BH" - P-value adjustment method (Benjamini-Hochberg FDR)
#   prevalence_filter: 0.002 - Minimum relative abundance threshold for filtering
#
# Required Packages:
#   igraph, psych, Hmisc, vegan, dplyr, reshape2
# ============================================================================

# Load required packages
library(igraph)
library(psych) 
library(Hmisc)
library(vegan) 
library(dplyr) 
library(reshape2) 

# Read OTU abundance table
# Format: rows = OTUs/ASVs, columns = samples, first column = row names
otu <- read.delim("network_RH_3h_1.txt", row.names = 1, check.names = FALSE)

# Filter genus-level data
# Retain OTUs/ASVs with relative abundance > 0.002 in at least 1 sample
otu <- otu[rowSums(otu > 0.002) >= 1, ]

# Calculate correlation matrix (R values and p-values)
cor_result <- rcorr(t(otu), type = "spearman")
occor.r = cor_result$r
occor.p = cor_result$P

# Adjust p-values using Benjamini-Hochberg FDR method
occor.p_adjusted <- p.adjust(occor.p, method = "BH")

# Filter significant correlations
# Set non-significant correlations to 0 based on p-value and correlation thresholds
occor.r[occor.p >= 0.05 | abs(occor.r) <= 0.75] = 0
diag(occor.r) <- 0 
occor.r[upper.tri(occor.r)] <- 0 

# Convert correlation matrix to long format
df = melt(as.matrix(occor.r)) 

# Construct Gephi-compatible edge file
df$Var1 = as.character(df$Var1) # Convert to character type 
df$Var2 = as.character(df$Var2)
df1 = subset(df, !df$Var1 == df$Var2) # Remove self-correlations 
colnames(df1) = c("Source", "Target", "Weight") 
df1 = subset(df1, !df1$Weight == 0) # Remove edges with zero weight

# Add absolute value of Weight column
df1$Weight_Abs <- abs(df1$Weight)

# Add Weight sign column (1 = positive, -1 = negative, 0 = zero)
df1$Sign <- ifelse(df1$Weight > 0, 1, ifelse(df1$Weight < 0, -1, 0))

# Export edge file
write.csv(df1, "221edge_RH_3h.csv", row.names = FALSE) 
# Read exported CSV file for column renaming
df <- read.csv("221edge_RH_3h.csv")

# Rename columns to standardized names
# Column order: Source, Target, original Weight, absolute Weight, correlation sign
names(df) <- c("Source", "Target",  "coll","Weight", "cor")

# Save updated edge file
write.csv(df, "221edge_RH_3h.csv", row.names = FALSE)

# Create node list from unique edge endpoints
df2 = data.frame(id = unique(c(df1$Source, df1$Target)))

# Read node attribute data (taxonomy or other metadata)
node_attributes <- read.delim("分类.txt", row.names = 1, check.names = FALSE)

# Merge node attributes with node list
# Filter and reorder attributes to match node list
node_attributes <- node_attributes[as.character(df2$id),]  # 按照节点列表的顺序进行筛选
df2 <- cbind(df2, node_attributes)  # 将节点属性合并到节点列表中

# Export final node file with attributes
write.csv(df2, "221node_with_attributes_RH_3h.csv", row.names = FALSE)

