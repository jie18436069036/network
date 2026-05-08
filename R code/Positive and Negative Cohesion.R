# ============================================================================
# Positive and Negative Cohesion
# Author: Yujie Wang
# Date: 2026-05-08
# 
# Description:
#   This script calculates positive and negative cohesion and connectedness 
#   metrics for microbial community analysis, following the framework by 
#   Herren & McMahon (2018). Cohesion quantifies the degree to which taxa 
#   within a community co-occur more or less often than expected by chance.
#   The script uses permutation-based null models to generate expected 
#   correlation matrices.
#
# Parameters:
#   input_file: "RH_3h.csv" - OTU/ASV abundance table (rows = samples, columns = taxa)
#   output_connectedness: "connectedness.csv" - Connectedness values per taxon
#   output_cohesion: "cohesion_RH_3h.csv" - Cohesion values per sample
#   pers.cutoff: 0.10 - Persistence cutoff (minimum fraction of samples a taxon must occupy)
#   iter: 200 - Number of permutations for null model (recommend >= 200)
#   tax.shuffle: TRUE - Use taxon shuffle (TRUE) or row shuffle (FALSE) for null model
#   use.custom.cors: FALSE - Whether to use a custom correlation matrix instead of calculating from data
#
# Output:
#   List containing four elements:
#     - Negative Connectedness: Mean negative observed-expected correlations per taxon
#     - Positive Connectedness: Mean positive observed-expected correlations per taxon
#     - Negative Cohesion: Sample-level negative cohesion scores
#     - Positive Cohesion: Sample-level positive cohesion scores
#
# Reference:
#   Herren, C.M. & McMahon, K.D. (2018). Cohesion: a method for quantifying 
#   the connectivity of microbial communities. The ISME Journal, 12, 793-802.
# ============================================================================

# Clear all variables from workspace
rm(list = ls())

# Count the number of zeros in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)}
# Calculate the mean of negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
  }
# Calculate the mean of positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
  }

# Set persistence cutoff (minimum presence proportion) for retaining taxa in analysis
pers.cutoff <- 0.10
# Set number of iterations for null model (recommend >= 200)
iter <- 200
# Choose taxon shuffle (tax.shuffle = TRUE) or row shuffle (tax.shuffle = FALSE)
tax.shuffle <- T

#Choose whether to use a custom correlation matrix
# Note: Your correlation table must have the same number of taxa as the abundance table.
# The abundance table should not have empty (all-zero) taxon vectors.
# Even if you input a custom correlation table, the persistence cutoff will still be applied.
use.custom.cors <- F

# Read dataset (each row is a sample)
b <-
  read.csv("RH_3h.csv", header = T, row.names = 1)
# If using custom correlation matrix, read and check dimensions
if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("RH_3h.csv", header = T,
                             row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  print(dim(b)[2] == dim(custom.cor.mat)[2])
  }

# Reformat data, remove empty samples and taxa
c <- as.matrix(b)
c <- c[rowSums(c) > 0,
       colSums(c) > 0]

# Save total individual counts for original samples
rowsums.orig <- rowSums(c)

# Determine zero count threshold for taxa based on persistence cutoff
zero.cutoff <- ceiling(pers.cutoff
                       * dim(c)[1])

# Remove taxa below persistence threshold
d <- c[ , apply(c, 2, zero) <
          (dim(c)[1]-zero.cutoff) ]

# Remove samples with no individuals
d <- d[rowSums(d) > 0, ]

# If using custom correlation matrix, update correlation matrix
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) <
                                         (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
  }

# Create relative abundance matrix
rel.d <- d / rowsums.orig

# View proportion of community retained
hist(rowSums(rel.d))

# Calculate observed correlation matrix
cor.mat.true <- cor(rel.d)

# Save median correlations for each taxon
med.tax.cors <- vector()

# Calculate expected correlations from null model
# Skip null model if using custom correlation matrix
if(use.custom.cors == F) {
  if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){
# Save correlations from each permutation
      perm.cor.vec.mat <- vector()
for(i in 1:iter){
  # Create empty matrix with same dimensions as rel.d
  perm.rel.d <- matrix(numeric(0),
                       dim(rel.d)[1], dim(rel.d)[2])
  rownames(perm.rel.d) <-
    rownames(rel.d)
  colnames(perm.rel.d) <-
    colnames(rel.d)
  # Shuffle each taxon
  for(j in 1:dim(rel.d)[2]){
    perm.rel.d[, j ] <- sample(rel.d[
      ,j ])
    }
  # Keep focal column unchanged
  perm.rel.d[, which.taxon] <- rel.d[
    , which.taxon]
  # Calculate correlation matrix of permuted matrix
  cor.mat.null <- cor(perm.rel.d)
  # Save correlations for focal taxon
  perm.cor.vec.mat <-
    cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
}
      # Save median correlations
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1,
                                                median))
if(which.taxon %% 20 == 0){print(which.taxon)}
      }
    } else {
      for(which.taxon in 1:dim(rel.d)[2]){
  # Save correlations from each permutation
        perm.cor.vec.mat <- vector()
  for(i in 1:iter){
    # Copy matrix for abundance shuffle
    perm.rel.d <- rel.d
    # Shuffle each sample
    for(j in 1:dim(rel.d)[1]){
      which.replace <- which(rel.d[j, ]
                             > 0 )
      which.replace.nonfocal <-
        which.replace[!(which.replace %in% which.taxon)]
    # Permute taxon vector
      perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal])
      }
    # Calculate correlation matrix of permuted matrix
    cor.mat.null <- cor(perm.rel.d)
    # Save correlations for focal taxon
    perm.cor.vec.mat <-
      cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
  }
        # Save median correlations
        med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1,
                                                  median))
  if(which.taxon %% 20 == 0){print(which.taxon)}
        }
      }
  }
# Calculate observed minus expected correlations
if(use.custom.cors == T) {
  obs.exp.cors.mat <- custom.cor.mat.sub
  } else {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
    }
diag(obs.exp.cors.mat) <- 0

# Calculate connectedness (mean of positive/negative correlations)
connectedness.pos <-
  apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <-
  apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion (relative abundance multiplied by connectedness)
cohesion.pos <- rel.d %*%
  connectedness.pos
cohesion.neg <- rel.d %*%
  connectedness.neg
# Output results 
output <- list(connectedness.neg,
               connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness","Positive Connectedness", "Negative Cohesion","Positive Cohesion")
print(output)

# Extract matrices or vectors to combine
connectedness_data <- list(
  "Negative Connectedness" = output[["Negative
                                     Connectedness"]],
  "Positive Connectedness" = output[["Positive
                                     Connectedness"]]
  )
# Convert to data frame and ensure same number of rows
connectedness_df <-
  as.data.frame(do.call(cbind, connectedness_data))

# Set column names
colnames(connectedness_df) <-
  names(connectedness_data)

# Write to CSV file
write.csv(connectedness_df, file ="connectedness.csv", row.names = TRUE)

# Extract matrices or vectors to combine
cohesion_data <- list(
  "Negative Cohesion" = output[["Negative Cohesion"]],
  "Positive Cohesion" = output[["Positive Cohesion"]])

# Convert to data frame and ensure same number of rows
cohesion_df <-
  as.data.frame(do.call(cbind, cohesion_data))
# Set column names
colnames(cohesion_df) <-
names(cohesion_data)
# Write to CSV file
write.csv(cohesion_df, file ="cohesion_RH_3h.csv", row.names = TRUE)
