Microbial Network Analysis

**Overview**

This repository contains the R code and data associated with the analysis of microbial co-occurrence networks and community cohesion.
Co-occurrence network analysis was performed to investigate microbial community interactions. Network robustness was assessed using natural connectivity metrics, and community cohesion was quantified using permutation-based null models.

**File Descriptions**
R Scripts

1.`network_plot.R` includes code to construct co-occurrence networks using Spearman correlation and generate edge/node files for network visualization.
- Data inputs:
  - `network_RH_3h_1.txt`: OTU/ASV abundance table (taxa as rows, samples as columns)
  - `Taxonomy.txt`: Taxonomic classification file for node attributes
- Data outputs:
  - `221edge_RH_3h.csv`: Network edge list with correlation weights and signs
  - `221node_with_attributes_RH_3h.csv`: Node list with taxonomic attributes
---
2.`natural_connectivity.R` includes code to assess network robustness by calculating natural connectivity after sequential random node removal.
- Data inputs:
  - `227edge_RH_3h.csv`: Network edge list from network construction
- Data outputs:
  - `natural_connectivity_results_RH_3h.csv`: Robustness scores at 10%, 20%, 30%, 40%, and 50% node removal levels
---
3.`Positive and Negative Cohesion.R` includes code to calculate positive and negative cohesion and connectedness following Herren & McMahon (2018).
- Data inputs:
  - `RH_3h.csv`: Sample-by-taxon abundance table (samples as rows, taxa as columns)
- Data outputs:
  - `connectedness.csv`: Taxon-level positive and negative connectedness values
  - `cohesion_RH_3h.csv`: Sample-level positive and negative cohesion scores

**Computing Environment**
- `Software：R ；Version：4.5.2；Purpose：Statistical computing environment
- `Software：RStudio ；Version：2025.05.1+513 (optional)；Purpose：Integrated development environment

**R Package Dependencies**
- `Package：igraph ；Version：2.2.2 ；Description：Network analysis and visualization 
- `Package：psych ；Version：2.6.1 ；Description：Correlation and descriptive statistics 
- `Package：Hmisc ；Version：5.2.5 ；Description：Data analysis utilities 
- `Package：vegan ；Version：2.7.3 ；Description：Community ecology analysis 
- `Package：dplyr ；Version：1.2.0 ；Description：Data manipulation
- `Package：reshape2 ；Version：1.4.5 ；Description：Data reshaping
- `Package：stats ；Version： 4.5.2 ；Description：Base R statistical functions

**Verify Package Versions**
```r
R.version.string
packageVersion("igraph")
packageVersion("psych")
packageVersion("Hmisc")
packageVersion("vegan")
packageVersion("dplyr")
packageVersion("reshape2")
