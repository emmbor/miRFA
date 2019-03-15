if (!dir.exists("results")){
  dir.create("results")
}

p <- argparser::arg_parser(description = "Description of program")

p <- argparser::add_argument(parser = p, "input", "input file")
p <- argparser::add_argument(parser = p, "--n", "name of run", default = NULL)
p <- argparser::add_argument(parser = p, "--g", "minimum number of genes need for enrichment", 
                             default = 5)
p <- argparser::add_argument(parser = p, "--t", "microt-CDS, should be between 0 and 1", 
                             default = 0.7)

argv <- argparser::parse_args(p)

source("pipeline/mirna_search_functions.R")

#Load packages
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))

mirnas<-read.table(argv$input,header=F)

run_name <- argv$n
min_n_genes <- argv$g
microt_cutoff <- argv$t


if (is.na(run_name)){
  run_name <- sub(pattern = "\\..*$", replacement = "", x = basename(argv$input))
}


whole_pipeline <- apply(X = mirnas, MARGIN = 1, FUN = mirna_search,
                       run_name = run_name,
                       min_n_genes = min_n_genes,
                       microt_cutoff = microt_cutoff)