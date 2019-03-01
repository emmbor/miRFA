setwd("/Users/emmyborgmastars/Documents/Rstats/databases")

#Run the mirna search script
source("mirna_search.R")

#Import list of mirnas
mirnas<-read.csv("mirna_list.txt",header=F)

#loop searches with each mirna using test function, margin = 1 (row), 2 (col)
whole_pipeline<-apply(X = mirnas, MARGIN = 1, FUN = mirna_search)

setwd("/Users/emmyborgmastars/Documents/Rstats/databases/DATA.DIR/Significant_correlations")

#Add number of miRNAs to 'vector'
vector<-c(1:15)

number_genes<-c()

#Create summary table with number of target genes
for (i in vector) {
  number_genes<-rbind(number_genes,whole_pipeline[[i]]$number_of_genes)
}

#All miRNA combined significant correlations

pos_mRNA<-c()
neg_mRNA<-c()
neg_prot<-c()
pos_prot<-c()

for (i in vector) {
  pos_mRNA<-rbind(pos_mRNA,whole_pipeline[[i]]$pos_temp_sign_mRNA)
}

for (i in vector) {
  neg_mRNA<-rbind(neg_mRNA,whole_pipeline[[i]]$neg_temp_sign_mRNA)
}

neg_prot<-c()
for (i in vector) {
  pos_prot<-rbind(pos_prot,whole_pipeline[[i]]$pos_temp_sign_prot)
}

for (i in vector) {
  neg_prot<-rbind(neg_prot,whole_pipeline[[i]]$neg_temp_sign_prot)
}

#Write to tables
#Save number of genes
setwd("/Users/emmyborgmastars/Documents/Rstats/databases/DATA.DIR")
write.table(number_genes,"number_of_genes.txt",sep="\t")
write.table(pos_mRNA,"pos_mRNA.txt",sep="\t")
write.table(neg_mRNA,"neg_mRNA.txt",sep="\t")
write.table(pos_prot,"pos_prot.txt",sep="\t")
write.table(neg_prot,"neg_prot.txt",sep="\t")

###########################################################################180324
#Check the correlations
setwd("/Users/emmyborgmastars/Documents/Rstats/databases")
par(mfrow=c(2,2))
check_cor<-function(X){
  gene<-X[2]
  miRNA<-X[1]
  
  mRNA_con <- dbConnect(SQLite(),'miRNAmRNACor.sqlite')
  ex_miRNA<-dbGetQuery(mRNA_con, paste0('SELECT * FROM miRNA WHERE miRNA_name IS"',   miRNA, '\n"'))
  v_miRNA<<-as.numeric(ex_miRNA[1,])
  
  ex_mRNA<-dbGetQuery(mRNA_con, 
                      paste0('SELECT * FROM mRNA WHERE hgnc_symbol IS"',   gene, '"'))
  v_mRNA<-as.numeric(ex_mRNA[1,])
  
  PCC_value<-cor(v_miRNA,v_mRNA,method="pearson",use="na.or.complete")
  plot(v_miRNA,v_mRNA,xlab=paste0(miRNA),ylab=paste0(gene),main=paste0('PCC: ',round(PCC_value, digits = 3)),
       cex=0.5,pch=1)
  abline(lm(v_mRNA~v_miRNA), col="red")
  
  dbDisconnect(mRNA_con)
}

#Check correlations
strong_pos_mRNA<-subset(pos_mRNA,PCC>=0.7,select=c(miRNA,Gene,PCC,P_value,adj_P_p))
strong_neg_mRNA<-subset(neg_mRNA,PCC<=-0.7,select=c(miRNA,Gene,PCC,P_value,adj_P_p))
strong_pos_prot<-subset(pos_prot,PCC>=0.7,select=c(miRNA,Protein,PCC,P_value,adj_P_p))
strong_neg_prot<-subset(neg_prot,PCC<=-0.7,select=c(miRNA,Protein,PCC,P_value,adj_P_p))

#Extract strong correlations
apply(X=strong_pos_mRNA,MARGIN = 1,FUN=check_cor)
apply(X=strong_neg_mRNA,MARGIN = 1,FUN=check_cor)
apply(X=strong_pos_prot,MARGIN = 1,FUN=check_cor)
apply(X=strong_neg_prot,MARGIN = 1,FUN=check_cor)

#Intersection of correlations on both mRNA and protein levels

for (i in vector) {
  setwd("/Users/emmyborgmastars/Documents/Rstats/databases/DATA.DIR")
  neg_prot<-rbind(neg_prot,whole_pipeline[[i]]$neg_temp_sign_prot)
  neg_prots<-whole_pipeline[[i]]$neg_temp_sign_prot
  pos_prots<-whole_pipeline[[i]]$pos_temp_sign_prot
  neg_mRNAs<-whole_pipeline[[i]]$neg_temp_sign_mRNA
  pos_mRNAs<-whole_pipeline[[i]]$pos_temp_sign_mRNA
  
  prots<-union(neg_prots$Protein,pos_prots$Protein)
  mRNAs<-union(neg_mRNAs$Gene,pos_mRNAs$Gene)
  both_levels<-intersect(prots,mRNAs)
  write.table(both_levels,paste0(i,'_both_levels.txt'),sep="\t")
}


