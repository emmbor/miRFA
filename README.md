# miRFA
Functional analysis of miRNA including 1) miRNA target prediction, 2) functional enrichment 3) correlation analyses between miRNA-target ge$

Description of programs:

mirna_search.R
- Contains the whole pipeline including miRNA target prediction, functional enrichment and correlation analysis.

masterpipeline.R
- The pipeline for calling the mirna_search.R script for 15 miRNAs and some post-processing.

3p_5p_annotation.R
- Function for summarizing 3p/5p miRNA isoform quantification for one file.

merge_mirna.R
- Function for merging several miRNA isoform quantification tables into one table.

3p5ppipeline.R
- Script for calling 3p_5p_annotation.R for all pancreatic adenocarcinoma-patients.
- Also for calling merge_mirna.R script in order to merge expression data for all patients.

MIMA_tr.pl
- Perl script for translating MIMAT-ID into miRNA names.
