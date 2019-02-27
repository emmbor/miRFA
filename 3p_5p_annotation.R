
library(plyr)

#function
convert_to_3p5p<-function(X){
    file_ID<-X[1]
    file_name<-X[2]
    sample_ID<-X[7]
  address<-(paste0('/Users/emmyborgmastars/Documents/Pythonstats/miRNA_isoform/', file_ID))
setwd(address)

data<-read.csv(file_name, header=TRUE,sep="\t")


#Extrahera ut endast rpm och miRNA
df<-data.frame(data$reads_per_million_miRNA_mapped,data$miRNA_region)
dfr<-subset(df,data.reads_per_million_miRNA_mapped >= 1)

#Summarize all rpm values of same miRNA identity
ny<-ddply(dfr,.(data.miRNA_region),summarize,sum_rpm=sum(data.reads_per_million_miRNA_mapped))

#Create log2 values, with sample_ID as colname
ny[sample_ID]<-log2(ny$sum_rpm)
colnames(ny)<-c("mature,miRNA_region","sum_rpm",sample_ID)

setwd("/Users/emmyborgmastars/Documents/Rstats/miRNA_data")
#Write to csv file
write.table(ny,paste0(sample_ID,'.csv'),sep=",")
setwd("/Users/emmyborgmastars/Documents/Rstats")
return(nrow(ny))


#Create  df



}


