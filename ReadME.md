# Project Structure:

# Objectives

This project aims to perform mutational signature analysis on cancer data, with the purpose of either dicovering mutational signatures from exposure or predicting mutational signatures from exisiting exposure. 

# Visualisations: 
- Top Exposures
- Proportional Exposures
- Exposures split by tumor types
- Exposures split by signature

## Comprehensive Downstream Analysis
- UMAP
- Coloured by Exposure

- HeatMap
- Compare Signatures
- Coloured by Clustering
- Differential Analysis

# Data
Mapping of TCGA variants by tumor type and signature exposure


# Steps

- i. Download TCGA data using TCGAbiolinks
- ii. Import TCGA samples into musicatk
- iii. Build mutational tables based on SBS96, DBS78 and IND83 schemas. 
- iv. Predict Sample Exposures using LDA
- v. Embed using UMAP
- vi. Coloured by tumor types and signatures


1. Data Types used: 
- MAF
- VCF 
- data.table 

2. Packages used: 
- musicatk

