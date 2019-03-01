#Load packages
library(DBI)
library(RSQLite)
library(VennDiagram)
library(edgeR)
library(GO.db)
library(org.Hs.eg.db)
path="/Users/emmyborgmastars/Documents/Rstats/databases"
setwd(path)

#Create a data frame number of genes
number_of_genes<-data.frame(matrix(ncol=4,nrow=0))
colnames(number_of_genes)<-c("miRNA_ID","Tarbase","microT_CDS","TargetScan")

#Create data frame for PCC 
PCC_table<-data.frame(matrix(ncol=4,nrow=0))
colnames(PCC_table)<-c("miRNA","Gene","PCC","P_value")

correlate_mRNA<-function(gene, miRNA){
mRNA_con <- dbConnect(SQLite(),'miRNAmRNACor.sqlite')
ex_mRNA<-dbGetQuery(mRNA_con, 
                  paste0('SELECT * FROM mRNA WHERE hgnc_symbol IS"',   gene, '"'))
v_mRNA<-suppressWarnings(as.numeric(ex_mRNA[1,]))

PCC_value<-suppressWarnings(cor(v_miRNA,v_mRNA,method="pearson",use="na.or.complete"))

if (!is.na(PCC_value)){
pval<-suppressWarnings(cor.test(v_miRNA,v_mRNA,method="pearson")$p.value)

#Add number of genes information to data frame
PCC_table<-rbind(PCC_table,
                  data.frame(miRNA=miRNA,Gene=gene,PCC=PCC_value,P_value=pval))
}
return(PCC_table)
dbDisconnect(mRNA_con)
}

#Create data frame for protein
PCC_p_table<-data.frame(matrix(ncol=4,nrow=0))
colnames(PCC_p_table)<-c("miRNA","Protein","PCC","P_value")

correlate_prot<-function(prot,miRNA){
  prot_con <- dbConnect(SQLite(),'miRNAProteinCor.sqlite')
  ex_prot<-dbGetQuery(prot_con, 
                      paste0('SELECT * FROM Protein WHERE protein_name IS"',   prot, '"'))

  v_prot<-suppressWarnings(as.numeric(ex_prot[1,]))

  PCC_value<-suppressWarnings(cor(p_miRNA,v_prot,method="pearson",use="na.or.complete"))
  
  if (!is.na(PCC_value)){
  p_val_s<-suppressWarnings(cor.test(p_miRNA,v_prot,method="pearson")$p.value)
  
  #Add number of genes information to data frame
  PCC_p_table<-rbind(PCC_p_table,
                    data.frame(miRNA=miRNA,Protein=prot,PCC=PCC_value,P_value=p_val_s))
  }
  return(PCC_p_table)
  dbDisconnect(prot_con)
}

mirna_search<-function(miRNA){

  #Connect to mirna_database
  con <- dbConnect(SQLite(),'mirna_database.sqlite')
  
  ##Search query mirna

  #Tarbase
  tarbase<-dbGetQuery(con, paste0('SELECT * FROM Tarbase WHERE mirna IS"',   miRNA, '"'))
  tarbase_g<-unique(tarbase$geneName)
  n_tarbase<-length(tarbase_g)
  
  #microT_CDS
  microT_CDS<-dbGetQuery(con, paste0('SELECT * FROM microT_CDS WHERE mirna IS"',   miRNA, '"','AND miTG_score>=0.7'))
  microT_g<-unique(microT_CDS$gene_name)
  n_microT<-length(microT_g)
  
  ###TargetScan
  #Conserved miRNA family
  targetscan1<-dbGetQuery(con, paste0('SELECT * FROM TargetScan_conserved WHERE miRNA IS"',   miRNA, '"'))
  targetscan1_g<-unique(targetscan1$Gene.Symbol)
  
  #Conserved site
  targetscan2<-dbGetQuery(con, paste0('SELECT * FROM TargetScan_conserved_site WHERE miRNA IS"',   miRNA, '"'))
  targetscan2_g<-unique(targetscan2$Gene.Symbol)
  
  #Nonconserved site
  targetscan3<-dbGetQuery(con, paste0('SELECT * FROM TargetScan_nonconserved_site WHERE miRNA IS"',   miRNA, '"'))
  targetscan3_g<-unique(targetscan3$Gene.Symbol)
  
  dbDisconnect(con)
  
  targetscan_g<-union(targetscan1_g,union(targetscan2_g,targetscan3_g))
  n_targetscan<-length(targetscan_g)
  write.table(targetscan_g,paste0('Target_genes/TargetScan_',miRNA,'.txt'))
  write.table(tarbase_g,paste0('Target_genes/Tarbase_',miRNA,'.txt'))
  write.table(microT_g,paste0('Target_genes/microT_CDS_',miRNA,'.txt'))
  
  ####Create Venn diagram
  
  #Calculate overlap
  overlap12 <- calculate.overlap(
    x = list(
      "Tarbase" = tarbase_g,
      "microT.CDS" = microT_g
    )
  );
  
  overlap23 <- calculate.overlap(
    x = list(
      "microT.CDS" = microT_g,
      "TargetScan" = targetscan_g
    )
  );
  
  overlap13 <- calculate.overlap(
    x = list(
      "Tarbase" = tarbase_g,
      "TargetScan" = targetscan_g
    )
  );
  
  overlap123 <- calculate.overlap(
    x = list(
      "Tarbase" = tarbase_g,
      "microT.CDS" = microT_g,
      "TargetScan" = targetscan_g
    )
  );
  
  area1<-length(tarbase_g)
  area2<-length(microT_g)
  area3<-length(targetscan_g)
  n12<-length(overlap12$a3)
  n23<-length(overlap23$a3)
  n13<-length(overlap13$a3)
  n123<-length(overlap123$a5)
  
  #create venndiagram with three sets
  
  color=c("thistle","slategray2","gainsboro")
  size=1.5
  grid.newpage()
  g=draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, 
                     category = c("Tarbase","microT-\n\t\t\t\t\t\t\tCDS","TargetScan"),
                     fill=color,cat.cex=rep(size,3),cex=rep(size,7))
  
  tiff(width = 7, height = 5, file=paste0('Venndiagrams/',miRNA,'.tiff'), units="in", res=600)
  grid.arrange(gTree(children=g), top=paste0(miRNA))
  dev.off()
  
  
  #Add number of genes information to data frame
  number_of_genes<-rbind(number_of_genes,
                    data.frame(miRNA_ID=miRNA,Tarbase=n_tarbase,
                    microT_CDS=n_microT,TargetScan=n_targetscan))

  #Union of target genes 
  uni<-union(tarbase_g,union(microT_g,targetscan_g))
  
  ###Correlation analyses

  #proteins
  #Create table
  prot_con <- dbConnect(SQLite(),'miRNAProteinCor.sqlite')
  e_miRNA<-dbGetQuery(prot_con, paste0('SELECT * FROM miRNA WHERE mirna_name IS"',   miRNA, '\n"'))
  p_miRNA<<-suppressWarnings(as.numeric(e_miRNA[1,]))
  cor_prot_list<-lapply(X=uni,FUN=correlate_prot,miRNA=miRNA)
  dbDisconnect(prot_con)
  
  #Multiple testing
  prot_df<-do.call(rbind.data.frame, cor_prot_list)
  adj_P_s<-p.adjust(prot_df$P_value, method="BH")
  prot_df<-cbind(prot_df,adj_P_s)
  temp_sign_prot<-subset(prot_df, adj_P_s < 0.05, select=c(miRNA,Protein,PCC,P_value,adj_P_s))
  pos_temp_sign_prot<-subset(temp_sign_prot, PCC > 0, select=c(miRNA,Protein,PCC,P_value,adj_P_s))
  neg_temp_sign_prot<-subset(temp_sign_prot, PCC <= 0, select=c(miRNA,Protein,PCC,P_value,adj_P_s))
  proteins<-union(as.vector(pos_temp_sign_prot$Protein),as.vector(neg_temp_sign_prot$Protein))
  
  #mRNAs
  #Create table 
  mRNA_con <- dbConnect(SQLite(),'miRNAmRNACor.sqlite')
  ex_miRNA<-dbGetQuery(mRNA_con, paste0('SELECT * FROM miRNA WHERE miRNA_name IS"',   miRNA, '\n"'))
  v_miRNA<<-suppressWarnings(as.numeric(ex_miRNA[1,]))
  cor_mRNA_list<-lapply(X=uni,FUN=correlate_mRNA,miRNA=miRNA)
  dbDisconnect(mRNA_con)
  
  #Multiple testing
  mRNA_df<-do.call(rbind.data.frame, cor_mRNA_list)
  adj_P_p<-p.adjust(mRNA_df$P_value, method="BH")
  mRNA_df<-cbind(mRNA_df,adj_P_p)
  temp_sign_mRNA<-subset(mRNA_df, adj_P_p < 0.05, select=c(miRNA,Gene,PCC,P_value,adj_P_p))
  pos_temp_sign_mRNA<-subset(temp_sign_mRNA, PCC > 0, select=c(miRNA,Gene,PCC,P_value,adj_P_p))
  neg_temp_sign_mRNA<-subset(temp_sign_mRNA, PCC <= 0, select=c(miRNA,Gene,PCC,P_value,adj_P_p))
  mRNAs<-union(as.vector(pos_temp_sign_mRNA$Gene),as.vector(neg_temp_sign_mRNA$Gene))
  
  uni2<-union(mRNAs,proteins)
  
  if(length(uni2)>0) {

  ###Functional enrichment after evaluation - union of mRNA and protein
  deID<-as.numeric(unlist(mapIds(org.Hs.eg.db,keys=uni2,
                                 column="ENTREZID",keytype="SYMBOL",multiVals="list")))
  
  #GO
  go<-goana(deID, FDR = 0.05, trend = FALSE, species="Hs")
  sign_GO<-subset(go,P.DE<=0.05,select=c(Term,Ont,N,DE,P.DE))
  GO_terms<-subset(sign_GO,DE>=5,select=c(Term,Ont,N,DE,P.DE))
  attach(GO_terms)
  write.table(GO_terms[order(-DE),],paste0('Functional_enrichment/After_evaluation/GO_',miRNA,'.txt'),
              sep=";")
  detach(GO_terms)
  
  #KEGG
  kegg<-kegga(deID, FDR = 0.05, trend = FALSE, species="Hs")
  sign_kegg<-subset(kegg,P.DE<=0.05,select=c(Pathway,N,DE,P.DE))
  kegg_paths<-subset(sign_kegg,DE>=5,select=c(Pathway,N,DE,P.DE))
  attach(kegg_paths)
  write.table(kegg_paths[order(-DE),],paste0('Functional_enrichment/After_evaluation/KEGG_',miRNA,'.txt'),
              sep=";")
  detach(kegg_paths)
  }

  return(list("number_of_genes"=number_of_genes,"pos_temp_sign_mRNA"=pos_temp_sign_mRNA,
              "neg_temp_sign_mRNA"=neg_temp_sign_mRNA,"pos_temp_sign_prot"=pos_temp_sign_prot,
              "neg_temp_sign_prot"=neg_temp_sign_prot))
}

