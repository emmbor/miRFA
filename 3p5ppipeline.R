setwd("/Users/emmyborgmastars/Documents/Rstats")

library(plyr)

##This package does not seem to be necessary
library(stringr)

#Run function script
source("3p_5p_annotation.R")

#Import data
df<-read.csv("gdc_sample_sheet_2018-02-13-03-23-11.tsv.2018-02-13T14_23_27.088749.tsv",sep="\t",header=T)

#create data frame
mirna_expr<-data.frame(matrix(ncol=1,nrow=1))
colnames(mirna_expr)<-"MIMA_ID"

#Run loop over each data row
apply(X = df, MARGIN = 1, FUN = convert_to_3p5p)

####MERGE
setwd("/Users/emmyborgmastars/Documents/Rstats")
source("merge_mirna.R")
setwd("/Users/emmyborgmastars/Documents/Rstats/miRNA_data")

m<-data.frame(matrix(ncol=1,nrow=1))
colnames(m)<-"mature.miRNA_region"


apply(X = df, MARGIN = 1, FUN = merge_mirna)

write.csv(m,"all_mirnas.csv")

#Transpose data frame
df<-read.csv("all_mirnas.csv",header=TRUE,sep=";")
df<-df[,-1]
df<-df[,-1]

tr<-t(df)

write.csv(tr,"t_mirnas.csv")

