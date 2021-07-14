
library(tidyverse)
library(DESeq2)
library(gplots)
library(ggdendro)
# library(pvclust)
library(limma)
# library(EnhancedVolcano)


#Load counts data.Use new consistent gene(2xy/(x^2+y^2)^0.5)
counts<-read.delim("sig_counts.txt",stringsAsFactors = FALSE,sep=" ")

# counts_new<-counts

genelist<-read.delim("new_score_consistent_genelist.txt",stringsAsFactors = FALSE)
counts_new<-left_join(genelist,counts,by = "GeneID")


counts_mat<-counts_new[,c(1:11)]%>%
  column_to_rownames("GeneID")%>%
  as.matrix()

# #load old conssitent gene
# counts<-read.delim("sig_counts.txt",stringsAsFactors = FALSE,sep=" ")
# genelist<-read.delim("old_score_consistent_genelist.txt",stringsAsFactors = FALSE)
# counts_new<-left_join(genelist,counts,by = "GeneID")
# counts_mat<-counts_new%>%
#   column_to_rownames("GeneID")%>%
#   as.matrix()



sampleinfo<-read.delim("Merged_Microglia_44_sampleInfo.txt",stringsAsFactors = FALSE)
sampleinfo<-sampleinfo[,c(1,2)]


#heatmap: use cpm counts.

#DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts_3group, colData=sampleinfo_3group, design=~CellType)
dds$CellType <- relevel(dds$CellType, ref = "Monocytes")
dds <- DESeq(dds,test="LRT")

# vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

#PCA w/ our samples
sampleinfo_3group<-sampleinfo[c(1:10),]
counts_3group<-counts_mat[,1:10]
#3group
dds_3group<-DESeqDataSetFromMatrix(countData=counts_3group, colData=sampleinfo_3group, design=~CellType)
dds_3group$CellType <- relevel(dds_3group$CellType, ref = "Monocytes")
dds_3group<-DESeq(dds_3group)

pc<- prcomp(t(cpm_mat))
tiff(filename="pca_with_ref.tiff",width = 1024, height = 960,res=150)
qplot(x=PC1,y=PC2,data=as.data.frame(pc$x),col=g$CellType,xlim = c(-50,30),ylim=c(-50,30),size=I(3))+
  labs(colour = 'cell type')
dev.off()


# |FC| (fold changes) >=1.2, P <= 0.001, and average expression >= 5 
res<-as.data.frame(results(dds,name="CellType_iMGCs_vs_Monocytes"))%>%
  rownames_to_column("GeneID")
res_filter<-res%>%
  filter(abs(log2FoldChange)>1.2&padj<1e-3&baseMean>5)
up<-res_filter%>%
  filter(log2FoldChange>0)
down<-res_filter%>%
  filter(log2FoldChange<0)
de_genelist<-data.frame("GeneID"=res_filter[,c(1)])

write.table(de_genelist,"deseq2_genelist.txt",quote=FALSE,row.names = FALSE)
  
# create heatmap.  
cpm_counts<-read.delim("consistenst_cpm_new.txt",stringsAsFactors = FALSE)  
# cpm_counts<-left_join(de_genelist,cpm_counts,by="GeneID")
cpm_mat<-cpm_counts%>%
  column_to_rownames("GeneID")%>%
  as.matrix()



group<-sampleinfo[,c(2,3)]
refinfo=data.frame(SampleID=colnames(cpm_counts)[12:45])
refinfo$CellType=c(rep("Ref_hMG",times=22),rep("Ref_Mono",times=12))
g=rbind(group,refinfo)
  
tiff(filename = "cluster_iMG_fixed.tiff",width = 2880, height = 2880,res=400)
heatmap.2(cpm_mat, col=redgreen(75), labCol=g$Sample_name,density.info="none", trace="none",scale="column",cexRow = 0.3,cexCol = 0.8,keysize = 1,dendrogram = "column")
dev.off()


# png(filename = "plot1.png")
# hplot1<-heatmap.2(counts.mat, scale="row", labCol=hgroups$cell.line, trace="none", density.info="none")
# invisible(dev.off())
# png("plot2.png")
# hplot2<-heatmap.2(t(counts.mat), col=redgreen(75), labRow=hgroups$cell.line,density.info="none", trace="none",scale="column")
# invisible(dev.off())
# png("plot3.png")
# hplot3<-heatmap.2(t(counts.mat), labRow=hgroups$cell.line,density.info="none", trace="none",scale="column")
# invisible(dev.off())



 

##################################################

  
#allcounts, venn diagram
  
   
allcounts<-read.delim("Original_mRNA_all_counts.txt",stringsAsFactors = FALSE)[,c(-1,-3,-4)]
nalist<-is.na(allcounts[,1])
allcounts<-allcounts[!nalist,]
 

 
#for all counts, remove duplicated
allcounts<-allcounts[!duplicated(allcounts),]
allcounts<-allcounts[order(allcounts$GeneSymbol),]
allcounts<-remove_rownames(allcounts)

repeated<-which(duplicated(allcounts[1]) | duplicated(allcounts[nrow(allcounts[1]):1, ])[nrow(allcounts[1]):1])
repeated<-sort(repeated,decreasing = TRUE)
backup<-allcounts
allcounts<-backup
for (i in repeated){
  if (apply(allcounts[i,c(2:11)], 1, function(x) all(x==0))){
    allcounts<-allcounts[-i,]
    print(i)
  } else if (apply(allcounts[i-1,c(2:11)], 1, function(x) all(x==0)) & allcounts[i-1,1]==allcounts[i,1]){
    d=i-1
    allcounts<-allcounts[-d,]
    print(d)
  } else if (apply(allcounts[i+1,c(2:11)], 1, function(x) all(x==0))& allcounts[i+1,1]==allcounts[i,1]) {
    d=i+1
    allcounts<-allcounts[-d,]
    print(d)
  } else{
    print(allcounts[i,1])
    allcounts<-allcounts[-i,]
  }
}

 


   
#for all counts
groups_all<-sampleinfo[c(1:10),]
groups_all$cell.line<- factor(groups_all$cell.line, 
                              levels = c("iMac","iMG", "monocyte"))
modelMatrix_all <- model.matrix(~ CellType, data = groups_all)
allcounts<-remove_rownames(allcounts)
countsmat_all<-allcounts%>%
  column_to_rownames("GeneSymbol")%>%
  as.matrix()



   
ddsObj <- DESeqDataSetFromMatrix(countData = countsmat_all,
                                 colData = groups_all,
                                 design = modelMatrix_all)
ddsObj <- DESeq(ddsObj)

counts_forpca<-t(counts(ddsObj, normalized=TRUE))
counts_forpca<-counts_forpca[,apply(counts_forpca, 2, function(x) all(x > 0))]
  
pcall<- prcomp(counts_forpca,center=TRUE,scale.=TRUE)
tiff(filename="pca_3group_allcounts.tiff",width = 1024, height = 960,res=150)
qplot(x=PC1,y=PC2,data=as.data.frame(pcall$x),col=sampleinfo_3group$CellType,size=I(3),xlim = c(-200,200),ylim=c(-200,200))+
  labs(colour = 'cell type')
dev.off() 


normcounts<-counts(ddsObj, normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("GeneID")
monocyte<-normcounts[,c(1,2,3,4)]
iMac<-normcounts[,c(1,5,6,7)]
iMG<-normcounts[,c(1,8,9,10,11)]

monocyte$mean<-rowMeans(monocyte[,-1])
iMac$mean<-rowMeans(iMac[,-1])
iMG$mean<-rowMeans(iMG[,-1])

monocytegene<-monocyte[monocyte$mean>=10,1]
iMacgene<-iMac[iMac$mean>=10,1]
iMGgene<-iMG[iMG$mean>=10,1]


write.table(monocytegene,"monocyte_counts.txt",quote=FALSE,row.names = FALSE)
write.table(iMacgene,"iMac_counts.txt",quote=FALSE,row.names = FALSE)
write.table(iMGgene,"iMG_counts.txt",quote=FALSE,row.names = FALSE)
cell.log<-cbind(monocyte,iMac,iMG)

countvenn <- vennCounts(cell.log)
vennDiagram(countvenn,cex=1,circle.col=c("red3", "blue2", "green3"))
tiff(filename="venn_highquality.tiff",width = 4000, height = 3200,res=400)
vennDiagram(countvenn,cex=0.8,circle.col=c("red3", "blue2", "green3"))

 

volcano
 
#for all gene
lfcres<-lfcShrink(ddsObj,coef = 2,independentFiltering=FALSE)
 
 
lfc<-lfcres


# lfc<-lfc[!is.na(lfc[,6]),]
# lfc<-lfc%>%
#   rownames_to_column("gene")%>%
#   filter(padj<=0.01)%>%
#   column_to_rownames("gene")
tiff(filename="volcano.tiff",width = 4500, height = 3000,res=400)
#volcanolist: padj<1E-3 and abs(logFC)>3
EnhancedVolcano(lfc,
                lab = rownames(lfc),
                x = 'log2FoldChange',
                y = 'pvalue',xlim = c(-5,5),ylim=c(0,50),pCutoff = 1e-3,selectLab = volcanolist)
 


Rearrange ORA result table and make barplot for it:  
   
ORA<-read.delim("ORA_results.txt",stringsAsFactors = FALSE)[,c(1:3)]
qv<-ORA[c(1,2,3,4,5,7,8,9,16,17),c(2,3)]
qv$pathway<-str_extract(qv$pathway,"^.+?(\\s|$).*?(\\s|$).*?(\\s|$)")
qv[5,2]<-"Cytokine Signaling"
qv[7,2]<-"MHC class II antigen"
qv[8,2]<-"Steroid biosynthesis"
qv[10,2]<-"Pyruvate metabolism and TCA cycle"
qv<-qv%>%
  mutate("log10qv"=-1*log10(qv$q.value))

tiff(filename="pathway_bar.tiff",width = 4000, height = 2000,res=400)
ggplot(data=qv,aes(x=reorder(pathway,log10qv),y=log10qv))+
  geom_bar(stat="identity",width=0.4)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="Pathway",y="-log10(q-value)")+
  coord_flip()
 


 
qvout<-ORA[c(1,2,3,4,5,7,8,9,16,17),c(3,2)]
write.table(qvout,file="Pathway_top10.txt",quote=FALSE,row.names = FALSE)
 


top pathway clustering  
Cholesterol:  
   
chel<-read.delim("chelosterol.txt",stringsAsFactors = FALSE,header=FALSE)
colnames(chel)<-"Geneid"
chel_counts<-left_join(chel,allcounts,by = c("Geneid" = "GeneSymbol"))
chel_counts<-chel_counts[,c(-2,-3,-4)]
chel_mat<-chel_counts%>%
  remove_rownames()%>%
  column_to_rownames("Geneid")%>%
  as.matrix()
groups_mac<-groups_all[c(4:10),]
tiff(filename="chel_heatmap.tiff",width = 4000, height = 3000,res=400)
heatmap.2(chel_mat, col=redgreen(75), labCol=groups_mac$cell.line,density.info="none", trace="none",scale="row",key.title = NA,Colv = FALSE,dendrogram ="none",keysize = 0.8,lhei=c(0.5,3),lwid=c(0.5,3),, key.par = list(cex=0.5), Rowv = TRUE)
 

immune
 
immu<-read.delim("immune.txt",stringsAsFactors = FALSE,header=FALSE)
colnames(immu)<-"Geneid"
immu_counts<-left_join(immu,allcounts,by = c("Geneid" = "GeneSymbol"))
immu_counts<-immu_counts[,c(-2,-3,-4)]
immu_mat<-immu_counts%>%
  remove_rownames()%>%
  column_to_rownames("Geneid")%>%
  as.matrix()
groups_mac<-groups_all[c(4:10),]
tiff(filename="immu_heatmap.tiff",width = 4000, height = 3000,res=400)
heatmap.2(immu_mat, col=redgreen(75), labCol=groups_mac$cell.line,density.info="none", trace="none",scale="row",key.title = NA,dendrogram ="none",keysize = 0.8,lhei=c(0.5,3),lwid=c(0.5,3),, key.par = list(cex=0.5), Rowv = TRUE)
 

cytokine
 
cyto<-read.delim("cytokine.txt",stringsAsFactors = FALSE,header=FALSE)
colnames(cyto)<-"Geneid"
cyto_counts<-left_join(cyto,allcounts,by = c("Geneid" = "GeneSymbol"))
cyto_counts<-cyto_counts[,c(-2,-3,-4)]
cyto_mat<-cyto_counts%>%
  remove_rownames()%>%
  column_to_rownames("Geneid")%>%
  as.matrix()
groups_mac<-groups_all[c(4:10),]
tiff(filename="cyto_heatmap.tiff",width = 4000, height = 3000,res=400)
heatmap.2(cyto_mat, col=redgreen(75), labCol=groups_mac$cell.line,density.info="none", trace="none",scale="row",key.title = NA,dendrogram ="none",keysize = 0.8,lhei=c(0.5,3),lwid=c(0.5,3),, key.par = list(cex=0.5), Rowv = TRUE)
 