# setwd("~/Nextcloud/Recherche/2_Plankton/Forams/2-Data_processing/3-Align_sequences")

# Process sequence alignment

Foram_Taxonomy <- readxl::read_xlsx("../../1-Data_raw/Data_Morard/Taxonomy/Table S5.xlsx", sheet = "Molecular Nomenclature")
head(Foram_Taxonomy)


## Spinose

Spinose_SSUalign_seq <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.95_pt0.95/SSU-align_Spinose/SSU-align_Spinose.eukarya.afa")
Spinose_SSUalign_seq <- Spinose_SSUalign_seq[order(Spinose_SSUalign_seq$seq.name),]
Spinose_SSUalign_seq
Spinose_SSUalign_seq_pf95_pt95 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.95_pt0.95/SSU-align_Spinose/SSU-align_Spinose.eukarya.mask.afa")
Spinose_SSUalign_seq_pf95_pt95 <- Spinose_SSUalign_seq_pf95_pt95[order(Spinose_SSUalign_seq_pf95_pt95$seq.name),]
Spinose_SSUalign_seq_pf95_pt95
Spinose_SSUalign_seq_pf5_pt5 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.5_pt0.5/SSU-align_Spinose/SSU-align_Spinose.eukarya.mask.afa")
Spinose_SSUalign_seq_pf5_pt5 <- Spinose_SSUalign_seq_pf5_pt5[order(Spinose_SSUalign_seq_pf5_pt5$seq.name),]
Spinose_SSUalign_seq_pf5_pt5

# Homogenize the typography for the unmask sequences
Spinose_SSUalign_seq$seq.text <- toupper(Spinose_SSUalign_seq$seq.text)
Spinose_SSUalign_seq$seq.text <- gsub("\\.", "-", Spinose_SSUalign_seq$seq.text)

# DNA instead of RNA
Spinose_SSUalign_seq$seq.text <- gsub("U", "T", Spinose_SSUalign_seq$seq.text)
Spinose_SSUalign_seq_pf95_pt95$seq.text <- gsub("U", "T", Spinose_SSUalign_seq_pf95_pt95$seq.text)
Spinose_SSUalign_seq_pf5_pt5$seq.text <- gsub("U", "T", Spinose_SSUalign_seq_pf5_pt5$seq.text)

# Match the nomenclature of modern lineages (MOTU2-level)
Foram_Taxonomy_Spinose <- Foram_Taxonomy[Foram_Taxonomy$Clade=="Spinose",]
Foram_Taxonomy_Spinose <- Foram_Taxonomy_Spinose[order(Foram_Taxonomy_Spinose$Taxonomic_Path),]
Foram_Taxonomy_Spinose

Spinose_SSUalign_seq$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Spinose_SSUalign_seq$seq.name)), "_", Foram_Taxonomy_Spinose$MOTU_lvl_1, Foram_Taxonomy_Spinose$MOTU_lvl_2)
Spinose_SSUalign_seq_pf95_pt95$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Spinose_SSUalign_seq_pf95_pt95$seq.name)), "_", Foram_Taxonomy_Spinose$MOTU_lvl_1, Foram_Taxonomy_Spinose$MOTU_lvl_2)
Spinose_SSUalign_seq_pf5_pt5$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Spinose_SSUalign_seq_pf5_pt5$seq.name)), "_", Foram_Taxonomy_Spinose$MOTU_lvl_1, Foram_Taxonomy_Spinose$MOTU_lvl_2)

# Keep only the first sequence for identical MOTU2 species
Spinose_SSUalign_seq_unique <- Spinose_SSUalign_seq[match(unique(Spinose_SSUalign_seq$seq.name), table=Spinose_SSUalign_seq$seq.name), ]
Spinose_SSUalign_seq_unique
Spinose_SSUalign_seq_pf95_pt95_unique <- Spinose_SSUalign_seq_pf95_pt95[match(unique(Spinose_SSUalign_seq_pf95_pt95$seq.name), table=Spinose_SSUalign_seq_pf95_pt95$seq.name), ]
Spinose_SSUalign_seq_pf95_pt95_unique
Spinose_SSUalign_seq_pf5_pt5_unique <- Spinose_SSUalign_seq_pf5_pt5[match(unique(Spinose_SSUalign_seq_pf5_pt5$seq.name), table=Spinose_SSUalign_seq_pf5_pt5$seq.name), ]
Spinose_SSUalign_seq_pf5_pt5_unique

phylotools::dat2fasta(Spinose_SSUalign_seq_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Spinose_SSU-align_alignment_species.fasta")
phylotools::dat2fasta(Spinose_SSUalign_seq_pf95_pt95_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Spinose_SSU-align_alignment_species_pf95_pt95.fasta")
phylotools::dat2fasta(Spinose_SSUalign_seq_pf5_pt5_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Spinose_SSU-align_alignment_species_pf5_pt5.fasta")


## Non-spinose

NonSpinose_SSUalign_seq <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.95_pt0.95/SSU-align_NonSpinose/SSU-align_NonSpinose.eukarya.afa")
NonSpinose_SSUalign_seq <- NonSpinose_SSUalign_seq[order(NonSpinose_SSUalign_seq$seq.name),]
NonSpinose_SSUalign_seq
NonSpinose_SSUalign_seq_pf95_pt95 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.5_pt0.5/SSU-align_NonSpinose/SSU-align_NonSpinose.eukarya.mask.afa")
NonSpinose_SSUalign_seq_pf95_pt95 <- NonSpinose_SSUalign_seq_pf95_pt95[order(NonSpinose_SSUalign_seq_pf95_pt95$seq.name),]
NonSpinose_SSUalign_seq_pf95_pt95
NonSpinose_SSUalign_seq_pf5_pt5 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.5_pt0.5/SSU-align_NonSpinose/SSU-align_NonSpinose.eukarya.mask.afa")
NonSpinose_SSUalign_seq_pf5_pt5 <- NonSpinose_SSUalign_seq_pf5_pt5[order(NonSpinose_SSUalign_seq_pf5_pt5$seq.name),]
NonSpinose_SSUalign_seq_pf5_pt5

# Homogenize the typography for the unmask sequences
NonSpinose_SSUalign_seq$seq.text <- toupper(NonSpinose_SSUalign_seq$seq.text)
NonSpinose_SSUalign_seq$seq.text <- gsub("\\.", "-", NonSpinose_SSUalign_seq$seq.text)

# DNA instead of RNA
NonSpinose_SSUalign_seq$seq.text <- gsub("U", "T", NonSpinose_SSUalign_seq$seq.text)
NonSpinose_SSUalign_seq_pf95_pt95$seq.text <- gsub("U", "T", NonSpinose_SSUalign_seq_pf95_pt95$seq.text)
NonSpinose_SSUalign_seq_pf5_pt5$seq.text <- gsub("U", "T", NonSpinose_SSUalign_seq_pf5_pt5$seq.text)

# Match the nomenclature of modern lineages (MOTU2-level)
Foram_Taxonomy_NonSpinose <- Foram_Taxonomy[Foram_Taxonomy$Clade=="Non-spinose",]
Foram_Taxonomy_NonSpinose <- Foram_Taxonomy_NonSpinose[order(Foram_Taxonomy_NonSpinose$Taxonomic_Path),]
Foram_Taxonomy_NonSpinose

NonSpinose_SSUalign_seq$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq$seq.name)), "_", Foram_Taxonomy_NonSpinose$MOTU_lvl_1, Foram_Taxonomy_NonSpinose$MOTU_lvl_2)
cbind(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq$seq.name)), Foram_Taxonomy_NonSpinose)
NonSpinose_SSUalign_seq_pf95_pt95$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq_pf95_pt95$seq.name)), "_", Foram_Taxonomy_NonSpinose$MOTU_lvl_1, Foram_Taxonomy_NonSpinose$MOTU_lvl_2)
cbind(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq_pf95_pt95$seq.name)), Foram_Taxonomy_NonSpinose)
NonSpinose_SSUalign_seq_pf5_pt5$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq_pf5_pt5$seq.name)), "_", Foram_Taxonomy_NonSpinose$MOTU_lvl_1, Foram_Taxonomy_NonSpinose$MOTU_lvl_2)
cbind(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", NonSpinose_SSUalign_seq_pf5_pt5$seq.name)), Foram_Taxonomy_NonSpinose)

# Keep only the first sequence for identical MOTU2 species
NonSpinose_SSUalign_seq_unique <- NonSpinose_SSUalign_seq[match(unique(NonSpinose_SSUalign_seq$seq.name), table=NonSpinose_SSUalign_seq$seq.name),]
NonSpinose_SSUalign_seq_unique
NonSpinose_SSUalign_seq_pf95_pt95_unique <- NonSpinose_SSUalign_seq_pf95_pt95[match(unique(NonSpinose_SSUalign_seq_pf95_pt95$seq.name), table=NonSpinose_SSUalign_seq_pf95_pt95$seq.name),]
NonSpinose_SSUalign_seq_pf95_pt95_unique
NonSpinose_SSUalign_seq_pf5_pt5_unique <- NonSpinose_SSUalign_seq_pf5_pt5[match(unique(NonSpinose_SSUalign_seq_pf5_pt5$seq.name), table=NonSpinose_SSUalign_seq_pf5_pt5$seq.name),]
NonSpinose_SSUalign_seq_pf5_pt5_unique

phylotools::dat2fasta(NonSpinose_SSUalign_seq_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/NonSpinose_SSU-align_alignment_species.fasta")
phylotools::dat2fasta(NonSpinose_SSUalign_seq_pf95_pt95_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/NonSpinose_SSU-align_alignment_species_pf95_pt95.fasta")
phylotools::dat2fasta(NonSpinose_SSUalign_seq_pf5_pt5_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/NonSpinose_SSU-align_alignment_species_pf5_pt5.fasta")


## Microperforate

Microperforate_SSUalign_seq <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.95_pt0.95/SSU-align_Microperforates/SSU-align_Microperforates.eukarya.afa")
Microperforate_SSUalign_seq <- Microperforate_SSUalign_seq[order(Microperforate_SSUalign_seq$seq.name),]
Microperforate_SSUalign_seq
Microperforate_SSUalign_seq_pf95_pt95 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.5_pt0.5/SSU-align_Microperforates/SSU-align_Microperforates.eukarya.mask.afa")
Microperforate_SSUalign_seq_pf95_pt95 <- Microperforate_SSUalign_seq_pf95_pt95[order(Microperforate_SSUalign_seq_pf95_pt95$seq.name),]
Microperforate_SSUalign_seq_pf95_pt95
Microperforate_SSUalign_seq_pf5_pt5 <- phylotools::read.fasta("../../3-Data_processed/Sequence_alignments/SSU-align/Mask_pf0.5_pt0.5/SSU-align_Microperforates/SSU-align_Microperforates.eukarya.mask.afa")
Microperforate_SSUalign_seq_pf5_pt5 <- Microperforate_SSUalign_seq_pf5_pt5[order(Microperforate_SSUalign_seq_pf5_pt5$seq.name),]
Microperforate_SSUalign_seq_pf5_pt5

# Remove species of independant evolutionary origin
Microperforate_SSUalign_seq <- Microperforate_SSUalign_seq[!grepl("Dentigloborotalia|Neogallitellia", Microperforate_SSUalign_seq$seq.name),]
Microperforate_SSUalign_seq_pf95_pt95 <- Microperforate_SSUalign_seq_pf95_pt95[!grepl("Dentigloborotalia|Neogallitellia", Microperforate_SSUalign_seq_pf95_pt95$seq.name),]
Microperforate_SSUalign_seq_pf5_pt5 <- Microperforate_SSUalign_seq_pf5_pt5[!grepl("Dentigloborotalia|Neogallitellia", Microperforate_SSUalign_seq_pf5_pt5$seq.name),]

# Homogenize the typography for the unmask sequences
Microperforate_SSUalign_seq$seq.text <- toupper(Microperforate_SSUalign_seq$seq.text)
Microperforate_SSUalign_seq$seq.text <- gsub("\\.", "-", Microperforate_SSUalign_seq$seq.text)

# DNA instead of RNA
Microperforate_SSUalign_seq$seq.text <- gsub("U", "T", Microperforate_SSUalign_seq$seq.text)
Microperforate_SSUalign_seq_pf95_pt95$seq.text <- gsub("U", "T", Microperforate_SSUalign_seq_pf95_pt95$seq.text)
Microperforate_SSUalign_seq_pf5_pt5$seq.text <- gsub("U", "T", Microperforate_SSUalign_seq_pf5_pt5$seq.text)

# Match the nomenclature of modern lineages (MOTU2-level)
Foram_Taxonomy_Microperforate <- Foram_Taxonomy[Foram_Taxonomy$Clade=="Microperforates",]
Foram_Taxonomy_Microperforate <- Foram_Taxonomy_Microperforate[order(Foram_Taxonomy_Microperforate$Taxonomic_Path),]
Foram_Taxonomy_Microperforate <- Foram_Taxonomy_Microperforate[!grepl("Dentigloborotalia|Neogallitellia", Foram_Taxonomy_Microperforate$Genus),]
Foram_Taxonomy_Microperforate

Microperforate_SSUalign_seq$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Microperforate_SSUalign_seq$seq.name)), "_", Foram_Taxonomy_Microperforate$MOTU_lvl_1, Foram_Taxonomy_Microperforate$MOTU_lvl_2)
Microperforate_SSUalign_seq_pf95_pt95$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Microperforate_SSUalign_seq_pf95_pt95$seq.name)), "_", Foram_Taxonomy_Microperforate$MOTU_lvl_1, Foram_Taxonomy_Microperforate$MOTU_lvl_2)
Microperforate_SSUalign_seq_pf5_pt5$seq.name <- paste0(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Microperforate_SSUalign_seq_pf5_pt5$seq.name)), "_", Foram_Taxonomy_Microperforate$MOTU_lvl_1, Foram_Taxonomy_Microperforate$MOTU_lvl_2)
cbind(gsub("\\|", "_", gsub("Representative_[0-9]*\\||Basetype_[0-9]*\\||basegroup_[0-9]*\\|", "", Microperforate_SSUalign_seq$seq.name)), Foram_Taxonomy_Microperforate)

# Keep only the first sequence for identical MOTU2 species
Microperforate_SSUalign_seq_unique <- Microperforate_SSUalign_seq[match(unique(Microperforate_SSUalign_seq$seq.name), table=Microperforate_SSUalign_seq$seq.name),]
Microperforate_SSUalign_seq_unique
Microperforate_SSUalign_seq_pf95_pt95_unique <- Microperforate_SSUalign_seq_pf95_pt95[match(unique(Microperforate_SSUalign_seq_pf95_pt95$seq.name), table=Microperforate_SSUalign_seq_pf95_pt95$seq.name),]
Microperforate_SSUalign_seq_pf95_pt95_unique
Microperforate_SSUalign_seq_pf5_pt5_unique <- Microperforate_SSUalign_seq_pf5_pt5[match(unique(Microperforate_SSUalign_seq_pf5_pt5$seq.name), table=Microperforate_SSUalign_seq_pf5_pt5$seq.name),]
Microperforate_SSUalign_seq_pf5_pt5_unique

phylotools::dat2fasta(Microperforate_SSUalign_seq_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Microperforate_SSU-align_alignment_species.fasta")
phylotools::dat2fasta(Microperforate_SSUalign_seq_pf95_pt95_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Microperforate_SSU-align_alignment_species_pf95_pt95.fasta")
phylotools::dat2fasta(Microperforate_SSUalign_seq_pf5_pt5_unique, outfile = "../../3-Data_processed/Sequence_alignments/SSU-align/Microperforate_SSU-align_alignment_species_pf5_pt5.fasta")
