setwd("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/Plots")

library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggsignif)
library(plyr)
library(ggpubr)
library(grid) 
library("pheatmap")
library("RColorBrewer")

####################################LOAD INPUT FILES

#Ampliconic gene copy number for 170 individuals (9 rows x 170 columns). The estimates are based read depth in from all sites on Y with mappability 1(First column of AmpliCoNE output).
AmpliconCN<-read.table("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/GTEx_CopyNumber_Amplicon.txt", sep="\t", header=T, stringsAsFactors = FALSE,row.names=1)
#colnames are GTEX-ID
colnames(AmpliconCN)<-substr(colnames(AmpliconCN),6,12)
rownames(AmpliconCN)[which(rownames(AmpliconCN)=="BPY")]="BPY2"
#XDG gene copy number for 170 individuals which are used in PCA analysis for removing outlier samples(12 rows x 170 columns)
XdgCN<-read.table("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/GTEx_CopyNumber_XDG.txt", sep="\t", header=T, stringsAsFactors = FALSE,row.names=1)
#colnames are GTEX-ID
colnames(XdgCN)<-substr(colnames(XdgCN),6,12)
rownames(XdgCN)[which(rownames(XdgCN)=="BPY")]="BPY2"
#rearranging column to make them similar
id_ordXDG<-match(colnames(AmpliconCN),colnames(XdgCN))
XdgCN<-XdgCN[,id_ordXDG]

#Family level gene expression values for 149 individuals post filtering(9 rows x 149 columns). DEseq2 normalized read counts.
GeneExpression<-read.table("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/results/GTEx_normalized_DESeq2_Amp.txt", sep="\t", header=T, stringsAsFactors = FALSE)
#colnames are GTEX-ID
colnames(GeneExpression)<-substr(colnames(GeneExpression),10,16)
rownames(GeneExpression)[which(rownames(GeneExpression)=="BPY")]="BPY2"
#Y Haplogroup information for each individual. Column:GTEX-ID, SRR-ID, HaplogroupInfo, HaplogroupInfo, Haplogroup. The last three column are parsed from Yhaplo out put but we only use the last column information for out analysis. 
Haplogroup<-read.table("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output/GTEX_Yhaplogroups_YchrBAM.tab", sep="\t", stringsAsFactors=FALSE,header=FALSE)
Yhapl<-cbind(substr(Haplogroup[,1],6,12),substr(Haplogroup[,5],1,1),substr(Haplogroup[,5],1,3))
colnames(Yhapl)<-c("Sample","Haplogroup","SubHaplogroup")

#Metadata for each individual, its a parsed file (first 5 column of the Phenotype.GRU file). Column : SUBJID, SEX, AGE, RACE, ETHNCTY.
Age<-read.table("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/results/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU_PARSED.txt", sep="\t", skip=10,stringsAsFactors=FALSE,header=TRUE,fill=TRUE)
metadata<-cbind(substr(Age$SUBJID,6,12),Age[,c(2:5)])
colnames(metadata)[1]<-c("Sample")

##############################FILTERING and GGPLOT friendly datastructure
Amp_Xdg<-rbind(AmpliconCN,XdgCN)
#Below steps were used to perform PCA analysis
# pca<-princomp(t(Amp_Xdg), cor=TRUE, scores=TRUE)
# biplot(pca)
# pdf("PCA_Amp_XDG_CN_allSamples_plot.pdf")
# plot(pca$scores[,1:2],type = "p",lwd=1, pch=19, main="PCA of XDG + Amp copy number")
# text(pca$scores[,1:2], labels = colnames(Amp_Xdg), pos = 4)
# dev.off()
#We removed only those file which have less than 1 copy number

outliers_samples<-c("NPJ8","11DXY","12BJ1")  #Samples whose ampliconic count is less than 1 or multiple XDG have copy numbers close to 0.5
outliers<-match(outliers_samples,colnames(Amp_Xdg))
Amp_Xdg_filtered<-Amp_Xdg[,-(outliers)]
AmpliconCN_filtered<-AmpliconCN[,-(outliers)]

#Creating index to combine all the input files.
Sample_Intersect<-intersect(colnames(GeneExpression),colnames(AmpliconCN_filtered))
idGE_orderSI<-match(Sample_Intersect,colnames(GeneExpression))
idCN_orderSI<-match(Sample_Intersect,colnames(AmpliconCN_filtered))
idYH_orderSI<-match(Sample_Intersect,Yhapl[,1])
idMD_orderSI<-na.omit(match(Sample_Intersect,metadata[,1]))

#Reording the columns in a standard order
AGE_filtered<-GeneExpression[,idGE_orderSI]
ACN_filtered<-AmpliconCN_filtered[,idCN_orderSI]
YH_filtered<-as.data.frame(Yhapl[idYH_orderSI,])
MD_filtered<-as.data.frame(metadata[as.vector(idMD_orderSI),])

#Creating ggplot friendly data structures 
GE_melt<-melt(t(AGE_filtered))  
CN_melt<-melt(t(ACN_filtered))
Yh_melt<-melt(YH_filtered)
##Summary of Copy Number and Gene Expression
if( all.equal(GE_melt$Var1,CN_melt$Var1) & all.equal(GE_melt$Var2,CN_melt$Var2)){ 
	Summary_AG<-cbind(GE_melt,CN_melt$value)
	colnames(Summary_AG)<-c("Sample","Gene","Expression","CopyNumber")
}

##Summary of all the inputs
if(all.equal(rep(Yh_melt[,1],9),GE_melt[,1])){
	Summary_AG_Yh<-cbind(Summary_AG,rep(Yh_melt[,2],9),rep(Yh_melt[,3],9))
	colnames(Summary_AG_Yh)<-c("Sample","Gene","Expression","CopyNumber","Haplogroup", "SubHaplogroup")
}

##Obtaining family specific data in the form of list
geneFamilies <- split( Summary_AG , f = Summary_AG$Gene )


#Colors used: 
colorPalette <- c("#E69F00", "#0078D7", "#009E73", "#D55E00", "#F0E442", "#56B4E9", "#CC79A7", "#9999CC", "#63a884")
# To use for fills, add
#scale_fill_manual(values=colorPalette)
# To use for line and point colors, add
#scale_colour_manual(values=colorPalette)


#########################################################Supplementary

medianAmpliconCN=log(apply(AmpliconCN_filtered,1,median))
varAmpliconCN=log(apply(AmpliconCN_filtered,1,var))
CNfit = lm(medianAmpliconCN ~ varAmpliconCN)
summary(CNfit)$r.squared
summaryAmpliconCN=cbind(rownames(AmpliconCN_filtered),medianAmpliconCN,varAmpliconCN)
colnames(summaryAmpliconCN)=c("Gene","Median","Variance")
ACN_M_sum=as.data.frame(summaryAmpliconCN,stringsAsFactors = FALSE)
ACN_M_sum$Median=as.numeric(ACN_M_sum$Median)
ACN_M_sum$Variance=as.numeric(ACN_M_sum$Variance)
pCNV<- ggplot(ACN_M_sum, aes(y=Variance, x=Median)) + geom_point(size=2,aes(color=Gene)) +
  geom_smooth(method=lm,se=FALSE)  +
  geom_point(size=4,aes(color=Gene))+
  #geom_text(size=6,aes(label=Gene,color=Gene),hjust=0, vjust=0, fontface="italic") +
  xlim(0, 3.8)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=12),
    axis.title.y = element_text(face="bold", color="black", size=12),
    plot.title = element_text(face="bold", color = "black", size=12),
    axis.text.x = element_text(face = "bold")) +
  #labs(title= "Copy number of Ampliconic Gene families") +
  labs(x = "log(Median copy number)") +
  labs(y = "log(Variance in copy number)")+
  annotate("text",x=3,y=-Inf,size=6,vjust=-1, hjust=-0.1,label=paste("R2=",round(summary(CNfit)$r.squared,digits = 2),by=""))+ scale_fill_manual(values=colorPalette)+ scale_colour_manual(values=colorPalette)
  ##theme(legend.position="none")

pdf("S1_GeneFamily_CopyNumber_MedianVariance_N167_plot.pdf")
pCNV
dev.off()

medianGeneExpression=log(apply(AGE_filtered,1,median))
varGeneExpression=log(apply(AGE_filtered,1,var))
GEfit = lm(medianGeneExpression ~ varGeneExpression)
summary(GEfit)$r.squared
summaryGeneExpression=cbind(rownames(AGE_filtered),medianGeneExpression,varGeneExpression)
colnames(summaryGeneExpression)=c("Gene","Median","Variance")
Exp_D_sum=as.data.frame(summaryGeneExpression,stringsAsFactors = FALSE)
Exp_D_sum$Median=as.numeric(Exp_D_sum$Median)
Exp_D_sum$Variance=as.numeric(Exp_D_sum$Variance)
pGEV <- ggplot(Exp_D_sum, aes(y=Variance, x=Median)) + geom_point(size=2,aes(color=Gene)) + geom_smooth(method=lm,se=FALSE)  +
  geom_point(size=4,aes(color=Gene))+
  #geom_text(size=6,aes(label=Gene,color=Gene),hjust=0, vjust=0, fontface="italic") +
  #xlim(3, 53)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(face="bold", color="black", size=12),
    axis.title.y = element_text(face="bold", color="black", size=12),
    plot.title = element_text(face="bold", color = "black", size=12),
    axis.text.x = element_text(face = "bold")) +
  #labs(title= "Expression of Ampliconic Gene families (VST)") +
  labs(x = "log(Median gene expression)") +
  labs(y = "log(Variance in gene expression)") +
  annotate("text",x=6,y=-Inf,size=6,vjust=-1, hjust=-0.1,label=paste("R2=",round(summary(GEfit)$r.squared,digits = 2),by=""))+ scale_fill_manual(values=colorPalette)+ scale_colour_manual(values=colorPalette)
  ##theme(legend.position="none")

pdf("S2_GeneFamily_Expression_MedianVariance_plot.pdf")
pGEV
dev.off()

###############GENE EXPRESS vs COPY NUM WITHIN FAMILY
GeneLevel_all_CNvsGE<-function(genelevel){
  r=cor.test(genelevel$CopyNumber,genelevel$Expression,method = "spearman")
  print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
  pGL<- ggplot(genelevel, aes(CopyNumber, Expression)) + 
    geom_point(aes(color=Gene))+
    geom_smooth(data=genelevel, mapping=aes(x = CopyNumber, y = Expression), col="gray", method=lm,se=FALSE) +
    labs(x = "Copy Number", y = "Expression", color = "Haplogroup\n")+
    theme_bw() + 
    theme(                              
      axis.title.x = element_text(face="bold",color="black"),
      axis.title.y = element_text(face="bold",color="black"),
      plot.title = element_text(face="bold.italic", color = "black"),
      legend.text=element_text(face="bold"),legend.title=element_text(face="bold"),legend.position="none") +
    labs(title= genelevel$Gene[1], hjust = 0) +
    labs(x = "Copy Number") +
    labs(y = "Gene Expression") +
    annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))#+
  #geom_text(size=2,aes(label=data_GH$Sample),hjust=0, vjust=0)
  return(pGL)
}

GeneLevel_CNvsGE<-function(genelevel,haplogroup){
  id=match(genelevel[,1],YH_filtered$Sample)
  #add haplogroup to gene level info
  data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
  data_GH<-data_GF[which(data_GF$SubHaplogroup==haplogroup),]
  r=cor.test(data_GH$CopyNumber,data_GH$Expression,method = "spearman")
  print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
  pGL<- ggplot(data_GH, aes(CopyNumber, Expression)) + 
    geom_point(aes(color=Gene))+
    geom_smooth(data=data_GH, mapping=aes(x = CopyNumber, y = Expression), col="gray", method=lm,se=FALSE) +
    labs(x = "Copy Number", y = "Expression", color = "Haplogroup\n")+
    theme_bw() + 
    theme(                              
      axis.title.x = element_text(face="bold",color="black"),
      axis.title.y = element_text(face="bold",color="black"),
      plot.title = element_text(face="bold.italic", color = "black"),
      legend.text=element_text(face="bold"),legend.title=element_text(face="bold"),legend.position="none") +
    labs(title= paste(data_GH$Gene[1],data_GH$SubHaplogroup[1],sep=" - "), hjust = 0) +
    labs(x = "Copy Number") +
    labs(y = "Gene Expression") +
    annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))#+
  #geom_text(size=2,aes(label=data_GH$Sample),hjust=0, vjust=0)
  return(pGL)
}
genelevel_plotsAll=llply(geneFamilies,GeneLevel_all_CNvsGE)
genelevel_plotsR1b=llply(geneFamilies,function(x){GeneLevel_CNvsGE(x,haplogroup="R1b")})
genelevel_plotsE1b=llply(geneFamilies,function(x){GeneLevel_CNvsGE(x,haplogroup="E1b")})
genelevel_plotsI1=llply(geneFamilies,function(x){GeneLevel_CNvsGE(x,haplogroup="I1a")})

pdf("S3_eneLevel_CNvsGE_ALL_plot.pdf")
grid.arrange(grobs=genelevel_plotsAll, ncol=3)
dev.off()

pdf("S4_GeneLevel_CNvsGE_haplR1b_plot.pdf")
grid.arrange(grobs=genelevel_plotsR1b, ncol=3)
dev.off()

pdf("S5_GeneLevel_CNvsGE_haplE1b_plot.pdf")
grid.arrange(grobs=genelevel_plotsE1b, ncol=3)
dev.off()

pdf("S6new_GeneLevel_CNvsGE_haplI1a_plot.pdf")
grid.arrange(grobs=genelevel_plotsI1, ncol=3)
dev.off()

###############GENE EXPRESSION VS AGE all samples

Age_CNvsGE<-function(genelevel){
	id=match(genelevel[,1],YH_filtered$Sample)
	genelevel_Yh<-cbind(genelevel,YH_filtered[id,c(2,3)])
	id1=match(MD_filtered$Sample,genelevel_Yh$Sample)
	genelevel_YhAge<-genelevel_Yh[id1,]
	AmpliconGECN_Age=cbind(as.matrix(genelevel_YhAge),as.vector(MD_filtered$AGE))
	colnames(AmpliconGECN_Age)[c(7)]=c("AGE")
	AmpliconGECN_Age<-as.data.frame(AmpliconGECN_Age,stringsAsFactors=F)
	AmpliconGECN_Age$AGE<-as.integer(AmpliconGECN_Age$AGE)
	AmpliconGECN_Age$Expression<-as.integer(AmpliconGECN_Age$Expression)
	AmpliconGECN_Age$CopyNumber<-as.integer(AmpliconGECN_Age$CopyNumber)
	r=cor.test(AmpliconGECN_Age$Expression,AmpliconGECN_Age$AGE,method = "spearman")
	print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
	pAGE1<-ggplot(AmpliconGECN_Age, aes(x=Expression,y=AGE)) +
		geom_point(aes(color=Gene))+ #aes(color=Sample)
		geom_smooth(method=lm,se=FALSE,col="gray")+
		annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		theme_bw() + 
		theme(                              
		  axis.title.x = element_text(face="bold", color="black"),
		  axis.title.y = element_text(face="bold", color="black"),
		  plot.title = element_text(face="bold.italic", color = "black")) +
		labs(title= paste(AmpliconGECN_Age$Gene[1])) +
		labs(x = "Gene Expression") +
		theme(legend.position="none")
	  #geom_text(size=2,aes(label=Sample),hjust=0, vjust=0)
	return (pAGE1)
}  

Age_CNvsGE_Haplogroup<-function(genelevel,haplogroup){
	id=match(genelevel[,1],YH_filtered$Sample)
	genelevel_Yh<-cbind(genelevel,YH_filtered[id,c(2,3)])
	id1=match(MD_filtered$Sample,genelevel_Yh$Sample)
	genelevel_YhAge<-genelevel_Yh[id1,]
	AmpliconGECN_Age=cbind(as.matrix(genelevel_YhAge),as.vector(MD_filtered$AGE))
	colnames(AmpliconGECN_Age)[c(7)]=c("AGE")
	AmpliconGECN_Age<-as.data.frame(AmpliconGECN_Age,stringsAsFactors=F)
	AmpliconGECN_Age$AGE<-as.integer(AmpliconGECN_Age$AGE)
	AmpliconGECN_Age$Expression<-as.integer(AmpliconGECN_Age$Expression)
	AmpliconGECN_Age$CopyNumber<-as.integer(AmpliconGECN_Age$CopyNumber)
	AmpliconGECN_Age<-AmpliconGECN_Age[(AmpliconGECN_Age$SubHaplogroup %in% haplogroup),]
	r=cor.test(AmpliconGECN_Age$Expression,AmpliconGECN_Age$AGE,method = "spearman")
	print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
	pAGE1<-ggplot(AmpliconGECN_Age, aes(x=Expression,y=AGE)) +
		geom_point(aes(color=Gene))+ #aes(color=Sample)
		geom_smooth(method=lm,se=FALSE,col="gray")+
		annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		theme_bw() + 
		theme(                              
		  axis.title.x = element_text(face="bold", color="black"),
		  axis.title.y = element_text(face="bold", color="black"),
		  plot.title = element_text(face="bold.italic", color = "black")) +
		labs(title= paste(AmpliconGECN_Age$Gene[1],AmpliconGECN_Age$SubHaplogroup[1],sep="-"),hjust = 0) +
		labs(x = "Gene Expression")+
		theme(legend.position="none")
	  #geom_text(size=2,aes(label=Sample),hjust=0, vjust=0)
	return (pAGE1)
}


familylevelGEvsAge_plots=llply(geneFamilies,Age_CNvsGE)
pdf("S6_GeneFamily_GeneExpression_AGE_ALL_plot.pdf")
grid.arrange(grobs=familylevelGEvsAge_plots, ncol=3)
dev.off()

familylevelGEvsAgeR1b_plots=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup(x,haplogroup="R1b")})
pdf("S7_GeneFamily_GeneExpression_AGE_R1b_plot.pdf")
grid.arrange(grobs=familylevelGEvsAgeR1b_plots, ncol=3)
dev.off()

familylevelGEvsAgeE1b_plots=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup(x,haplogroup="E1b")})
pdf("S8_GeneFamily_GeneExpression_AGE_E1b_plot.pdf")
grid.arrange(grobs=familylevelGEvsAgeE1b_plots, ncol=3)
dev.off()

familylevelGEvsAgeI1a_plots=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup(x,haplogroup="I1a")})
pdf("S9new_GeneFamily_GeneExpression_AGE_I1a_plot.pdf")
grid.arrange(grobs=familylevelGEvsAgeI1a_plots, ncol=3)
dev.off()