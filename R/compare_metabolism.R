library(plyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(lattice)
library(Vennerable)
#library(VennDiagram)
library(colorspace)
#source("http://bioconductor.org/biocLite.R"); biocLite(c("graph", "RBGL"))
#install.packages("reshape")
#install.packages("Vennerable", repos="http://R-Forge.R-project.org")

compare_meta <- read.csv('/home/stephan/test/results/analysis/compare_meta_3.csv',sep="\t",header = TRUE)
df <- data.frame(compare_meta[,2:5] , row.names=compare_meta[,2])

#ggplot(data=df,aes(df["933-7","JBM10","JBC07","JBNZ41"]))+ geom_col(data=df)
#ggplot(data=df)+ geom_col(data=df,aes(y = df[,1], df[,2]))
###?names(df)
#help(package = "Vennerable")
#?barplot

vek <- data.matrix(df[,2:4])
rownames(vek) <- rownames(df)
barplot(vek,main ="Gene count of a metabolic pathway assortment",
        xlab="counts",ylab = "strains",legend.text = FALSE,horiz = FALSE, 
        beside = TRUE)
# transponierte matrix
#tvek <- t(vek)
# col=topo.colors(3)
par(mar=c(14,8,1,7))
x <- barplot(t(vek),axisnames=FALSE,legend.text = TRUE, cex.names = 0.3,
            args.legend = list(x = 'right',xpd=TRUE,bty='n', inset=c(-0.15,0),y.intersp= 3,pt.cex = 0.01),beside = TRUE,horiz = FALSE,ylab="counts",col=c("red","royalblue","yellowgreen"))#c("#FFFF00","#A0522D","#0000FF"),"#FFFF00","#A0522D","#0000FF","#FFFF00","green4","yellowgreen"

#x <- barplot(t(vek)[2:6,],axisnames=FALSE,legend(x = 'right',legend= colnames(vek),xpd=TRUE,bty='n', inset=c(-0.15,0)), cex.names = 1.0,
#        beside = TRUE,horiz = FALSE,ylab="counts", col=topo.colors(5))
x <- x[1,]
labs <- paste(rownames(vek))
#x <- barplot(table(mtcars$cyl), xaxt="n")
#labs <- paste(names(table(mtcars$cyl)), "cylinders")
text(cex=1, x=x+2.25, y=-10, labs, xpd=TRUE, srt=55, adj= 1)






# Ven-diagram ko_ids
#ko_dict <- vector(mode="character", length=0)
#for (species in c(933-7,JBM10,JBC07,JBNZ41,DS)){
#  now <- read.csv(paste('/home/stephan/test/results/analysis/',species,'_ko-list.csv',sep = ""),sep="\t",header = FALSE,strip.white =TRUE)
#  ko_id <- t(now[1])
#  ko_dict <- cbind(ko_dict,c(ko_id))
#}
#var933 <-t(read.csv('/home/stephan/test/results/analysis/933-7_ko-list.csv',sep="\t",header = FALSE)[1])
#JBM10 <-t(read.csv('/home/stephan/test/results/analysis/blast/JBM10_similar_genes_all.count',sep="\t",header = FALSE)[1])
#JBC07 <-t(read.csv('/home/stephan/test/results/analysis/blast/JBC07_similar_genes_all.count',sep="\t",header = FALSE)[1])
#JBNZ41 <-t(read.csv('/home/stephan/test/results/analysis/blast/JBNZ41_similar_genes_all.count',sep="\t",header = FALSE)[1])
JBM10 <-t(read.csv('/home/stephan/test/results/extracted_genes/JBM10_unique_genes.txt',sep="\t",header = FALSE)[1]) 
JBC07 <-t(read.csv('/home/stephan/test/results/extracted_genes/JBC07_unique_genes.txt',sep="\t",header = FALSE)[1])
JBNZ41 <-t(read.csv('/home/stephan/test/results/extracted_genes/JBNZ41_unique_genes.txt',sep="\t",header = FALSE)[1])
#DS <-t(read.csv('/home/stephan/test/results/analysis/DS_ko-list.csv',sep="\t",header = FALSE)[1])
Poteriospumella <- unique(c(JBNZ41,JBC07,JBM10))
#VennDiagrams(Vstem3, doWeights = FALSE,circle.col = c("red", "blue", "green3"))
#?Venn
#x <- list( JBM10_ = c(unique(JBM10)),JBC07_ =c(unique(JBC07)),JBNZ41_ = c(unique(JBNZ41))) #, var933_7 = c(var933),DS_ = c(DS)

x <- list( JBM10_ = c(JBM10),JBC07_ =c(JBC07),JBNZ41_ = c(JBNZ41)) #, var933_7 = c(var933),DS_ = c(DS)
Vstem <- Venn(x)
plot(Vstem, doWeights = TRUE)
plot(Vstem, doWeights = FALSE)
#plot(Vstem, doWeights = TRUE, show = list(FaceText = "signature", SetLabels = FALSE,  Faces = TRUE))


y <- list(JBM10_ = c(JBM10) ,JBC07_ =c(JBC07),JBNZ41_ = c(JBNZ41))
#y <- list(P.lacustris = c(Poteriospumella) ,P.malhamensis = c(DS),O.danica = c(var933))
Vstem3 <- Venn(y)
gp <- VennThemes(Vstem3,colourAlgorithm="sequential")
plot(Vstem3,gp=gp,doWeights = FALSE)
#plot(Vstem3, doWeights = FALSE)

?plot
#heatmap
heat_matrix <- matrix(data = NA,nrow=5,ncol = 5)
species <- list(unique(c(var933)),unique(c(JBM10)),unique(c(JBC07)),unique(c(JBNZ41)),unique(c(DS)))
rownames(heat_matrix) <- c("933-7","JBM10","JBC07","JBNZ41","DS")
colnames(heat_matrix) <- rownames(heat_matrix)
for (xrow in 1:5){a
  for (ycol in 1:(6-xrow)){
    min_gene_number <- min(
      length(species[[xrow]]),length(species[[ycol]])
      )
    percent <- length((intersect(species[[xrow]],species[[ycol]])))/min_gene_number
    heat_matrix[xrow,ycol]<- round(percent,digits=2)
  }
}

d2_matrix <- melt(heat_matrix)
windows(title="Similarity analysis",16,11)
p <- ggplot(d2_matrix, aes(Var1,Var2)) +
  geom_tile(aes(fill = value), colour = "white")+
  scale_fill_gradientn(colours=heat.colors(5), name="No. clusters")+ ###oder3
  geom_text(aes(Var1,Var2, label = value), color = "black", size=5)+
  theme(legend.title = element_text(colour="black", size=14, face="bold"))+
  theme(legend.text=element_text(size=14))+
  theme(axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=14, face="bold"))
p

# assembly, gene analysis - statistic
statistic_1 <-read.csv('/home/stephan/Documents/statistic_assembly_gene_analysis.csv',sep="\t",header = TRUE)
stat_mat <- data.matrix(statistic_1)
stat_df <- data.frame(statistic_1, row.names=1)
rownames(stat_mat) <- rownames(stat_df)
par(mar=c(8,5,5,5))
barplot(stat_mat[,2:6],axisnames=FALSE,legend.text = TRUE, cex.names = 0.3, log = "y", names.arg = colnames(stat_df),
        args.legend = list(x = 'right','top',xpd=TRUE,bty='n', inset=c(-0.15,0),y.intersp= 1,pt.cex = 0.01),beside = TRUE,horiz = FALSE,ylab="counts", col=c("bisque1","bisque3","bisque4","slategray2","darkblue"))

##test
