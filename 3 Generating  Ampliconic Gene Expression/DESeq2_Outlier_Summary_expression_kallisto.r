/galaxy/home/rxv923/.Renv/versions/3.3.0/bin/R

library(tximport)
library("DESeq2")
tar2gen<-read.table("/galaxy/home/rxv923/refs/Hg/hg38/RNA/hg38_refFlat_10252016_GenetoID.txt", header=F, sep="\t", stringsAsFactors=F)
names(tar2gen)=c("target_id","gene")

base_dir <- "/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/results/Kallisto_Results"
samples<-list.dirs(base_dir,recursive = FALSE,full.names=FALSE)
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id,"abundance.tsv"))

#Using tximport to read kallisto output files
gene_level=tximport(kal_dirs, type = "kallisto",countsFromAbundance = "no", tx2gene = tar2gen, reader = read.delim, geneIdCol="gene",txIdCol="target_id" )

#Obtaining the gene level counts 
gtexData=gene_level$counts
gtexData<-as.matrix(gtexData)
gtexData=apply(gtexData, 1:2, round)
head(gtexData)

#DESeq2 steps to set the design of the data. No replicates and all are paired end datasets 
condition <- colnames(gtexData)
type <- rep("paired-end",length(condition))

colData=data.frame(cbind(condition,type),stringsAsFactors = FALSE)
row.names(colData)=colnames(gtexData)
head(colData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gtexData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

d<-counts(dds, normalized=TRUE)
r<-assay(rld)
v<-assay(vsd)

save(dds,d,v,r, file = "DESeq_gtex_normalized.RData")
write.table(d,file="GTEx_Counts_Testis_171samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(r,file="GTEx_RLDCounts_Testis_171samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(v,file="GTEx_VSTCounts_Testis_171samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#########################################################################################################
library(ggplot2)
library(gridExtra)
library("DESeq2")

load("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/results/DESeq_gtex_normalized.RData")
head(d)
head(r)
head(v)

sampleID<-substr(colnames(d),10,15)
pca_v<-prcomp(t(log(v)))
plot(pca_v, type = "l")

pdf("PCA_Testis_GE_allSamples_plot.pdf")
plot(pca_v$x[,1:2], type = "p",lwd=4, main= "PCA of testis samples normalized gene expression")
#text(pca_v$x[,1:2], labels = sampleID, pos = 4)
dev.off()
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#Generate Heat map of the count matrix
library("pheatmap")
library("RColorBrewer")
#Select the top 1000 genes based on mean values of the row
select <- order(rowMeans(v),decreasing=TRUE)[1:1000]

#Plot a heatmap of the top 1000 genes 
pdf("Heatmap_Testis_top1000_GE_allSamples_plot.pdf")
pheatmap(v[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE)
dev.off()


#apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists <- dist(t(v))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sampleID
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Heatmap_Testis_bwsample_distance_GE_allSamples_plot.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,show_rownames=FALSE)
dev.off()


		 
		 
Correlation_matrix<-cor(v)
#Here first row is kallisto_111CU which is not an outlier
#outliers<-which(Correlation_matrix[1,]<0.95)
outliers_samples<-c("kallisto_117YX","kallisto_1192W","kallisto_11DXY","kallisto_11UD2","kallisto_12BJ1","kallisto_131XF","kallisto_13G51","kallisto_13OW8","kallisto_13QJ3","kallisto_145MF","kallisto_NPJ8", "kallisto_P4PQ","kallisto_PLZ6","kallisto_QDVN","kallisto_S4Z8","kallisto_SIU7","kallisto_WY7C","kallisto_X4XX","kallisto_ZLFU","kallisto_ZY6K","kallisto_ZDTT")
outliers<-match(outliers_samples,colnames(v))

v_o<-v[,-(outliers)]
sampleID_o<-substr(colnames(v_o),10,15)
pca_vo<-prcomp(t(log(v_o)))
plot(pca_vo, type = "l")
pdf("PCA_Testis_GE_Filtered_plot.pdf")
plot(pca_vo$x[,1:2], type = "p",lwd=4, main= "PCA of testis samples normalized gene expression")
text(pca_vo$x[,1:2], labels = sampleID, pos = 4)
dev.off()

select_o <- order(rowMeans(v_o),decreasing=TRUE)[1:1000]
pdf("Heatmap_Testis_top1000_GE_Filtered_plot.pdf")
pheatmap(v_o[select_o,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE)
dev.off()
sampleDists_o <- dist(t(v_o))
sampleDistMatrix_o <- as.matrix(sampleDists_o)
rownames(sampleDistMatrix_o) <- sampleID_o
colnames(sampleDistMatrix_o) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Heatmap_Testis_bwsample_distance_GE_Filtered_plot.pdf")
pheatmap(sampleDistMatrix_o,
         clustering_distance_rows=sampleDists_o,
         clustering_distance_cols=sampleDists_o,
         col=colors,show_rownames=FALSE)
dev.off()	 

############################################################################################################################
#Post filtered dataset
outliers_samples<-c("kallisto_117YX","kallisto_1192W","kallisto_11DXY","kallisto_11UD2","kallisto_12BJ1","kallisto_131XF","kallisto_13G51","kallisto_13OW8","kallisto_13QJ3","kallisto_145MF","kallisto_NPJ8", "kallisto_P4PQ","kallisto_PLZ6","kallisto_QDVN","kallisto_S4Z8","kallisto_SIU7","kallisto_WY7C","kallisto_X4XX","kallisto_ZLFU","kallisto_ZY6K","kallisto_ZDTT")
outliers<-match(outliers_samples,colnames(v))
d_final<-d[,-(outliers)]
v_final<-v[,-(outliers)]
genelist=c("SRY","AMELY","DDX3Y","EIF1AY","NLGN4Y","PRKY","KDM5D","TBL1Y","TMSB4Y","USP9Y","UTY","ZFY")
ampliconList=c("BPY2B","BPY2","BPY2C","CDY1","CDY1B","CDY2A","CDY2B","DAZ1","DAZ2","DAZ3","DAZ4","HSFY1","HSFY2","PRY","PRY2","RBMY1A1","RBMY1J","RBMY1F","RBMY1E","RBMY1D","RBMY1B","TSPY10","TSPY1","TSPY3","TSPY8","TSPY4","TSPY2","VCY","VCY1B","XKRY","XKRY2")

xdg=NULL
for (gene in genelist){
  print(gene)
  xdg=rbind(xdg,d_final[which(row.names(d_final)==gene),])
}
xdg_normD=cbind(genelist,xdg)


amp=NULL
for (gene in ampliconList){
  print(gene)
  amp=rbind(amp,d_final[which(row.names(d_final)==gene),])
}
amp_normD=cbind(ampliconList,amp)


Amplicon_summary=rbind(
  BPY=apply(amp[1:3,],2,sum),
  CDY=apply(amp[4:7,],2,sum),
  DAZ=apply(amp[8:11,],2,sum),
  HSFY=apply(amp[12:13,],2,sum),
  PRY=apply(amp[14:15,],2,sum),
  RBMY=apply(amp[16:21,],2,sum),
  TSPY=apply(amp[22:27,],2,sum),
  VCY=apply(amp[28:29,],2,sum),
  XKRY=apply(amp[30:31,],2,sum)
)


write.table(xdg_normD,file="GTEx_normalized_DESeq2_XDG.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(Amplicon_summary,file="GTEx_normalized_DESeq2_Amp.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

