library(edgeR)
myData <- read.table("All_data_merge_FPKM_197.txt", head=T, sep="\t")
# names(myData)
# length(myData$gene_id)

my.mat<-myData%>%
  column_to_rownames("gene_id")%>%
  as.matrix()

#remove batch effect
my.mat<-removeBatchEffect(my.mat,batch=sampleinfo$Batch)


# transpose the dataset
myData.T <- as.data.frame(t(my.mat))


# make a duplicate
myData.T2 <- myData.T


# assign batch numbers
myData.T2$Batch<-sampleinfo$Batch


# start calculation for mds following the instruction
# at this website:
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/122-multidimensional-scaling-essentials-algorithms-and-r-code/#demo-data

library(magrittr)
library(dplyr)
library(ggpubr)
library(ggplot2)

# rearrange the data columns so that the sample name and batch are in the first two columns
myData.T2 <- subset(myData.T2, select=c(15437, 1:15436))
bckupt2=myData.T2

myData.T2<-bckupt2


rnum.keep=which(sampleinfo$Batch %in% c(1,4,7,8))
celltype.rm=which(grepl("CD44|iPS|NPC",sampleinfo$CellType))
rnum.keep=setdiff(rnum.keep,celltype.rm)
group.keep=sampleinfo[rnum.keep,]
myData.T2=myData.T2[rnum.keep,]


# # view the columns to verify
# myData.T2[1:10, 1:10]

# start calculation for mds
mds <- myData.T2 %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
head(mds)

# batch <- as.factor(MDS$Batch)
# MDS$CellType <- as.factor(MDS$CellType)
# ggplot(mds, aes(Dim.1, Dim.2))+ geom_point(aes(shape=as.factor(group.keep$CellType), color=as.factor(group.keep$Batch)))                           
#                          

# 
# # plot mds
# ggscatter(mds, x = "Dim.1", y = "Dim.2",
#           size = 1,
#           repel = TRUE)
# 
# 
# 
# 
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           # label = rownames(myData.T2),
#           size = 1,
#           repel = TRUE)
# 
# 
# 
# 
# 
# 
# myData$Batch <- as.factor(myData$Batch)
# myData.T2$Batch <- as.factor(myData.T2$Batch)
# # Compute MDS
# mds <- myData.T2 %>%
#   dist() %>%          
#   cmdscale() %>%
#   as_tibble()
# colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           # label = rownames(myData.T2),
#           size = 1,
#           repel = TRUE)
# # K-means clustering
# clust <- kmeans(mds, 3)$cluster %>%
#   as.factor()
# mds <- mds %>%
#   mutate(groups = clust)
# # Plot and color by groups
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           # label = rownames(myData.T2),
#           color = "groups",
#           palette = "jco",
#           size = 1, 
#           ellipse = TRUE,
#           ellipse.type = "convex",
#           repel = TRUE)
# ls()
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(myData.T2),
#           color = "groups",
#           palette = "jco",
#           size = 1, 
#           ellipse = TRUE,
#           ellipse.type = "convex",
#           repel = TRUE)
# ls()

# Compute MDS for myData, which does not have the batch numbers
# Compute MDS
bckup=myData.T

myData.T=bckup
myData.T=myData.T[rnum.keep,]

mds <- myData.T %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           # label = rownames(myData.T),
#           size = 1,
#           repel = TRUE)
# ls()
# head(mds)
# mds <- myData.T %>%
#   dist() %>%          
#   cmdscale()
# head(mds)
# plot(mds)
# mds <- as.data.frame(mds)
# mds$samples <- myData.T2$Batch
# head(mds)
names(mds[1]) <- "Dim_1"
names(mds[2]) <- "Dim_2"
# head(mds)
# names(mds)[2] <- "Dim_2"
# names(mds)[1] <- "Dim_1"
# head(mds)

# FMG_Muffat
# iMGL_Muffat
# MGL_Gosselin
# Mono_Gosselin
# iMGL_Grubman


#change legend
group.keep$label<-c(rep("FMG_Muffat",time=3),rep("iMGL_Muffat",time=9),rep("MGL_Gosselin",time=22),rep("Mono_Gosselin",time=12),rep("iMGL_Grubman",time=13),rep("Mono",time=3),rep("iMacs",time=3),rep("iMGCs",time=4))

p <- ggplot(mds, aes(Dim_2,Dim_1)) + geom_point(aes(shape=as.factor(group.keep$CellType), color=as.factor(group.keep$Batch)))+
  scale_color_manual(name="Batch",labels=c("Muffat","Gosselin","Grubman","JC"),
                     values=c("red","darkgreen","blue","purple"))+
  scale_shape_manual(name="CellType",labels=c("MGL","FMG","iMacs","iMGCs","monocytes","iMGL"),values=c(0:5))

tiff(filename="mds_change_legend.tiff",width=2000,height=2000,res=300)                     
p

#add label in the  plot
p <- ggplot(mds, aes(Dim_2,Dim_1)) + geom_point(aes(shape=as.factor(group.keep$CellType), color=as.factor(group.keep$Batch)))+
  geom_text(label=group.keep$label,nudge_x = 0.25, nudge_y = 0.25,check_overlap = T )+
  scale_color_manual(name="Batch",labels=c("Muffat","Gosselin","Grubman","JC"),
                     values=c("red","darkgreen","blue","purple"))+
  scale_shape_manual(name="CellType",labels=c("MGL","FMG","iMacs","iMGCs","monocytes","iMGL"),values=c(0:5))
tiff(filename="mds_add_lable_on_plot.tiff",width=3000,height=3000,res=300)   
p
dev.off()
