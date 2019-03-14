if (!dir.exists("results")){
  dir.create("results")
}

p <- argparser::arg_parser(description = "Description of program")

p <- argparser::add_argument(parser = p, "input", "input file")
p <- argparser::add_argument(parser = p, "--n", "name of run", default = NULL)

argv <- argparser::parse_args(p)

source("mirna_search_functions.R")

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

if (is.na(run_name)){
  run_name <- sub(pattern = "\\..*$", replacement = "", x = basename(argv$input))
}
print(getwd())
whole_pipeline <- apply(X = mirnas, MARGIN = 1, FUN = mirna_search, run_name = argv$n)