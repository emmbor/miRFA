#Function for merging

merge_mirna<-function(X){
  sample_ID<-X[7]
  
  #Import data
  x<-read.csv(paste0(sample_ID,'.csv'),header=T)
  x<-data.frame(x[1],x[3])
  m<<-merge(x,m,all=TRUE,by="mature.miRNA_region")
}