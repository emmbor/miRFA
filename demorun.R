# DEMO
#############################################################################
# Add your path
path = "/Users/emmyborgmastars/desktop/demo"

setwd(path)

mirnas<-read.table("demo_input.txt",header=F)

# Optional - add name of run


name_run <- function(mirnas, name_of_run = NULL){
  if (is.null(name_of_run)){
    name_of_run <- deparse(substitute(mirnas))
  }
  return(name_of_run)
}


runname<-name_run(mirnas,name_of_run="Run1")

#Run the mirna search script

source("mirna_search.R")

#loop searches with each mirna using test function, margin = 1 (row), 2 (col)
whole_pipeline<-apply(X = mirnas, MARGIN = 1, FUN = mirna_search)
#############################################################################








