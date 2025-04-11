# Overview
This project focuses on studying the diversification dynamics of planktonic foraminifera using various macroevolutionary models. Key goals include understanding cryptic speciation and leveraging both fossil and molecular data to reconstruct evolutionary patterns. The models used include Protracted Birth-Death (PBD), Occurrence Birth-Death Diffusion (OBDD), and GeoSSE.

Some files are too large for this repository, including the Triton database (which can be found at https://doi.org/10.6084/m9.figshare.c.5242154.v4), and others that can be regenerated with `2-Data_processing` scripts.

# Directory Structure and Description
## 1-Data_raw
Contains the raw datasets, including:

- `Aze_phylogeny`: Morphospecies phylogenetic dataset from Aze et al.
- `Cenozoic_depth_rasters`: Raster files providing depth information for the Cenozoic era to analyze environmental context.
- `Data_Morard`: Genetic and taxonomic data for extant species from Raphaël Morard.
- `Functional_data`: Dataset of functional traits.
- `TritonDB`: Occurrence records from the Triton database.

## 2-Data_processing
Scripts and intermediate files for processing raw data into a format suitable for downstream analyses:

- `1-Process_Triton_occurrences`: Scripts to clean, preprocess and subsample TritonDB occurrences for phylogenetic and biogeographic analyses.
- `2-Process_Aze_phylogeny`: Steps to process the phylogenies from Aze et al.
- `3-Align_sequences`: Sequence alignment workflows using structural information on the ribosomal RNA (SSU-align).
- `4-Assign_biogeography`: Assign biogeographical regions to species based on their occurrence records.

## 3-Data_processed
Processed datasets ready for analysis:

- `Biogeography`: Final biogeographical data with assigned regions.
- `Cenozoic_depth_rasters_plots`: Visualizations of depth rasters for the Cenozoic era, used to correlate with phylogenetic patterns.
- `Morphospecies_phylogenies`: Generated phylogenies for morphospecies in different lineages.
- `Sequence_alignments`: Aligned sequence data using SSU-align.
- `Triton_occurrences`: Cleaned and preprocessed occurrences from TritonDB.

## 4-Phylogenetic_reconstruction
Final outputs from phylogenetic reconstruction, including maximum clade credibility (MCC) trees and other dated phylogenies.

## 5-Analyses
Scripts and outputs for running different diversification models:

- `1-ProSSE`: Analysis using the Protracted State-Dependent Speciation and Extinction (ProSSE) model to study cryptic diversity.
- `2-GeoSSE`: GeoSSE analysis results, investigating diversification dependent on biogeographic ranges.
- `3-ESSE`: Extinction and Speciation State-dependent models (ESSE) results.
- `4-OBDD`: Occurrence Birth-Death Diffusion (OBDD) model outputs using fossil occurrences and modern molecular data.

# Replication Steps
1. Run Data Processing Scripts: Execute scripts in 2-Data_processing to clean and process the raw data. These scripts prepare occurrence data, perform sequence alignment with SSU-align, and assign biogeographical regions.
2. Build Phylogenies: Use Rev scripts in 4-Phylogenetic_reconstruction to perform phylogenetic inference.
3. Diversification Analysis: Run the diversification models (ProSSE, GeoSSE, ESSE, and OBDD) using the scripts in 5-Analyses. Adjust model parameters as needed based on the data processed and phylogenetic tree generated.

# Softwares
- R (with packages: diversitree, hisse, ape, phytools, dplyr, ggplot2, ggtree, paleoPhylo, divDyn, readxl, etc.)
- SSU-align for sequence alignment
- RevBayes for phylogenetic inference using OBDP
- Julia for running Tapestree.jl

# Contact Information
For any questions or issues replicating this pipeline, please reach out to:

Jérémy Andréoletti: jeremy.andreoletti@gmail.com
