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


####################FIGURES

#Colors used: 
colorPalette <- c("#E69F00", "#0078D7", "#009E73", "#D55E00", "#F0E442", "#56B4E9", "#CC79A7", "#9999CC", "#63a884")
# To use for fills, add
#scale_fill_manual(values=colorPalette)
# To use for line and point colors, add
#scale_colour_manual(values=colorPalette)


#1
###Function to invert the X axis labels for better readability
draw_colnames_90 <- function (coln, gaps, ...)
{
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,
        "bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
    return(res)
}
#Annotation of AZF regions
#AZFb includes seven ampliconic gene families (BPY2, CDY, DAZ, HSFY, PRY, RBMY, and XKRY) and
#AZFc which partially overlaps with AZFb, includes five such families (BPY2, CDY, DAZ, PRY, and RBMY)
#annotation_col = data.frame(AZF = factor(rep(c("AZFb_c", "AZFb", "AZFb_c","nonAZF", "AZFb"), c(3, 1, 2,2,1))))
#rownames(annotation_col)<-rownames(AGE_filtered)
#annotation_colors = list(AZF = c(AZFb_c = "#009E73", AZFb="#F0E442",nonAZF="white"))
#Loading the draw_colnames_90 function and plotting heatmap for copy number and gene expression
#Copy Number
assignInNamespace(x="draw_colnames", value="draw_colnames_90",ns=asNamespace("pheatmap"))
cnhm<-pheatmap(cor(t(AmpliconCN_filtered)),cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE,main="A",color = colorRampPalette(c("#0078D7","white", "#D55E00"))(70), treeheight_row=0,treeheight_col=0) #annotation_col = annotation_col,annotation_colors = annotation_colors
#Gene expression
assignInNamespace(x="draw_colnames", value="draw_colnames_90",ns=asNamespace("pheatmap"))
gehm<-pheatmap(cor(t(AGE_filtered)), cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE,main="B",treeheight_row=0,treeheight_col=0,color = colorRampPalette(c("#0078D7", "white", "#D55E00"))(70)) #, annotation_col = annotation_col,annotation_colors = annotation_colors

#Print
pdf("Fig1new_Correlation_GeneFamily_GEandCN_plot_minusAZFregion.pdf")
grid.arrange(grobs=list(cnhm[[4]],gehm[[4]]),nrow=2)
dev.off()


#2
#Convert the Copy number and gene expression to log scale
logAmpliconVCN=Summary_AG
logAmpliconVCN$CopyNumber=log(Summary_AG$CopyNumber+0.01) #Adding 0.01 to prevent 0 conversion
logAmpliconVCN$Expression=log(Summary_AG$Expression+1) #Adding 1 to prevent 0 conversion. 1 and not 0.01 because gene expression counts are large numbers and they can handle addition of 1.

#Liner regression model to caliculate the R-square value.
CNvsGEfit = lm(logAmpliconVCN$CopyNumber ~ logAmpliconVCN$Expression)
summary(CNvsGEfit)$r.squared

#Plot
pCVvGEwithLine<- ggplot(logAmpliconVCN, aes(x = CopyNumber, y = Expression) )+ 
  geom_boxplot(size=1,outlier.shape=NA,aes(color=factor(Gene)),position=position_dodge(width = 0.5)) +
  geom_point(size=0.7,aes(color=factor(Gene)))+
  geom_smooth(data=logAmpliconVCN, mapping=aes(x = CopyNumber, y = Expression), col="black",method=lm,se=FALSE) +
  labs(x = "Copy Number", y = "Expression", color = "Gene families\n",shape="Gene families\n")+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(face="bold",color="black", size=12),
    axis.title.y = element_text(face="bold",color="black", size=12),
    plot.title = element_text(face="bold", color = "black", size=12),
    legend.text=element_text(face="bold.italic"),legend.title=element_text(face="bold")) +
  #labs(title= "Ampliconic Copy Number vs Expression") +
  labs(x = "Ampliconic Gene Copy Number (log scale)") +
  labs(y = "log(Ampliconic Gene Expression)") + 
  ##Convert the log scale to natural scale in X axis for Copy number
  scale_x_continuous(breaks=c(0,0.6931472,1.098612,1.609438,2.079442,2.484907,2.995732,3.496508),labels=c("0" = "1","0.6931472" = "2", " 1.098612" = "3","1.609438"="5","2.079442"="8", "2.484907" = "12", "2.995732" = "20", "3.496508" = "33"))

#Coordinates where the family names should be printed. Customized postions without significance only for better view
x=c(1.09, 1.90, 1.55, 0.46, 0.50, 2.50, 3.25, 0.30, 0.35 )
y=c(4.75, 5.733, 7.88, 6.12, 3.50, 6.25, 8.5, 7.90, 1.80)
lab=unique(logAmpliconVCN$Gene)

#Print
pdf("Fig2_GeneFamily_Expression_CopyNumber_plot.pdf")
pCVvGEwithLine + annotate("text", x = x, y = y, label = lab,fontface = 'italic',size=6)+ theme(legend.position="none")+ scale_fill_manual(values=colorPalette)+ scale_colour_manual(values=colorPalette)
dev.off()






#3
#Function for group level copy number and gene expression across Y haplogroups. The margins are adjusted so CN and GE can be plotted below each other
GL_YhvsGE<-function(genelevel){
	#add haplogroup to gene level info
	id=match(genelevel[,1],YH_filtered$Sample)
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	#parse data with respect to the four major haplogroups 
	data_GF<-data_GF[(data_GF$Haplogroup %in% c("E","R","I","J")),]
	#Order in which the haplogroups will be printed on the X-axis
	data_GF$Haplogroup<- factor(data_GF$Haplogroup, levels = c("E","J","I","R"))
	c1<-ggplot(data_GF, aes(x=Haplogroup, y=Expression))+
		geom_violin(aes(color=Haplogroup))+
		geom_boxplot(aes(color=Haplogroup),outlier.shape=NA,width=0.1) +
		geom_hline(yintercept = mean(data_GF$Expression), linetype = 2)+
		stat_compare_means(method = "anova",aes(label = paste0("p =", ..p.format..)),label.y=max(data_GF$Expression)+0.1*max(data_GF$Expression))+
		#stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",hide.ns = TRUE)+
		theme_bw() + 
		theme(axis.title.x = element_blank(),axis.text.x = element_text(),plot.title = element_blank())+
		labs(title = data_GF$Gene[1])+
		labs(y = "Expression", color = "Haplogroup\n") +
		#margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
		theme(plot.margin=unit(c(1,1,1,1),unit = "pt"))+
		theme(legend.position="none")+ scale_fill_manual(values=colorPalette)+ scale_colour_manual(values=colorPalette)
	return(c1)
}


GL_YHvsCN<-function(genelevel){
	#add haplogroup to gene level info
	id=match(genelevel[,1],YH_filtered$Sample)
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	#parse data with respect to the four major haplogroups 
	data_GF<-data_GF[(data_GF$Haplogroup %in% c("E","R","I","J")),]
	#Order in which the haplogroups will be printed on the X-axis
	data_GF$Haplogroup<- factor(data_GF$Haplogroup, levels = c("E","J","I","R"))
	c1<-ggplot(data_GF, aes(x=Haplogroup, y=CopyNumber))+
		geom_violin(aes(color=Haplogroup))+
		geom_boxplot(aes(color=Haplogroup),outlier.shape=NA,width=0.1) +
		geom_hline(yintercept = mean(data_GF$CopyNumber), linetype = 2)+
		stat_compare_means(method = "anova",aes(label = paste0("p =", ..p.format..)),label.y=max(data_GF$CopyNumber)+0.1*max(data_GF$CopyNumber))+
		#stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",hide.ns = TRUE)+
		theme_bw() + 
		theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),plot.title = element_text(face="bold.italic", color = "black", size=12,hjust = 0.5))+
		#margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
		theme(plot.margin=unit(c(1,1,0,1),unit = "pt"))+
		labs(title= data_GF$Gene[1])+
		labs(y = "Copy Number", color = "Haplogroup\n") +
		theme(legend.position="none")+ scale_fill_manual(values=colorPalette)+ scale_colour_manual(values=colorPalette)
	return(c1)
}

#Calling the function on each gene family
familylevelExpressionvsHaplo_plots=llply(geneFamilies,GL_YhvsGE)
familylevelCopyNumbervsHaplo_plots=llply(geneFamilies,GL_YHvsCN)


#Arranging the subplots in order of 3 x 3 gene families
g1A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$BPY2))
g1B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$BPY2))
g2A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$CDY))
g2B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$CDY))
g3A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$DAZ))
g3B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$DAZ))
g4A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$HSFY))
g4B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$HSFY))
g5A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$PRY))
g5B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$PRY))
g6A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$RBMY+ ylim(5,17.5)))
g6B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$RBMY))
g7A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$TSPY+ ylim(20,85)))
g7B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$TSPY))
g8A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$VCY))
g8B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$VCY))
g9A=ggplot_gtable(ggplot_build(familylevelCopyNumbervsHaplo_plots$XKRY))
g9B=ggplot_gtable(ggplot_build(familylevelExpressionvsHaplo_plots$XKRY))

#Making sure all the plots have same width
maxWidth = grid::unit.pmax(g1A$widths[2:3], g1B$widths[2:3],g2A$widths[2:3], g2B$widths[2:3],g3A$widths[2:3], g3B$widths[2:3],g4A$widths[2:3],g4B$widths[2:3],g5A$widths[2:3],g5B$widths[2:3],g6A$widths[2:3],g6B$widths[2:3],g7A$widths[2:3],g7B$widths[2:3],g8A$widths[2:3],g8B$widths[2:3],g9A$widths[2:3],g9B$widths[2:3])
g1A$widths[2:3] <- as.list(maxWidth)
g1B$widths[2:3] <- as.list(maxWidth)
g2A$widths[2:3] <- as.list(maxWidth)
g2B$widths[2:3] <- as.list(maxWidth)
g3A$widths[2:3] <- as.list(maxWidth)
g3B$widths[2:3] <- as.list(maxWidth)
g4A$widths[2:3] <- as.list(maxWidth)
g4B$widths[2:3] <- as.list(maxWidth)
g5A$widths[2:3] <- as.list(maxWidth)
g5B$widths[2:3] <- as.list(maxWidth)
g6A$widths[2:3] <- as.list(maxWidth)
g6B$widths[2:3] <- as.list(maxWidth)
g7A$widths[2:3] <- as.list(maxWidth)
g7B$widths[2:3] <- as.list(maxWidth)
g8A$widths[2:3] <- as.list(maxWidth)
g8B$widths[2:3] <- as.list(maxWidth)
g9A$widths[2:3] <- as.list(maxWidth)
g9B$widths[2:3] <- as.list(maxWidth)

#Print
grid.newpage()
pdf('Fig3new_Haplogroup_GEandCN_2.pdf', width=6.5, height=9)
grid.arrange(arrangeGrob(g1A,g2A,g3A,g1B,g2B,g3B,g4A,g5A,g6A,g4B,g5B,g6B,g7A,g8A,g9A,g7B,g8B,g9B,ncol=3,heights=rep(0.163,6))    )
dev.off()

#4
load("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/results/DESeq_gtex_normalized.RData")

outliers_samples<-c("kallisto_117YX","kallisto_1192W","kallisto_11DXY","kallisto_11UD2","kallisto_12BJ1","kallisto_131XF","kallisto_13G51","kallisto_13OW8","kallisto_13QJ3","kallisto_145MF","kallisto_NPJ8", "kallisto_P4PQ","kallisto_PLZ6","kallisto_QDVN","kallisto_S4Z8","kallisto_SIU7","kallisto_WY7C","kallisto_X4XX","kallisto_ZLFU","kallisto_ZY6K","kallisto_ZDTT")
outliers<-match(outliers_samples,colnames(v))
d<-d[,-(outliers)]


ampliconList=c("BPY2B","BPY2","BPY2C","CDY1","CDY1B","CDY2A","CDY2B","DAZ1","DAZ2","DAZ3","DAZ4","HSFY1","HSFY2","PRY","PRY2","RBMY1A1","RBMY1J","RBMY1F","RBMY1E","RBMY1D","RBMY1B","TSPY10","TSPY1","TSPY3","TSPY8","TSPY4","TSPY2","VCY","VCY1B","XKRY","XKRY2","HSFX1","HSFX2","VCX","VCX2","VCX3A","VCX3B","DAZL","CDYL","CDYL2","TSPYL2","RBMX","RBMX2","RBMXL3")

amp=NULL
for (gene in ampliconList){
  print(gene)
  amp=rbind(amp,d[which(row.names(d)==gene),])
}
amp_normV=cbind(ampliconList,amp)
Amplicon_summary=rbind(
  BPY2_Family=apply(amp[1:3,],2,sum),
  CDY_Family=apply(amp[4:7,],2,sum),
  DAZ_Family=apply(amp[8:11,],2,sum),
  HSFY_Family=apply(amp[12:13,],2,sum),
  PRY_Family=apply(amp[14:15,],2,sum),
  RBMY_Family=apply(amp[16:21,],2,sum),
  TSPY_Family=apply(amp[22:27,],2,sum),
  VCY_Family=apply(amp[28:29,],2,sum),
  XKRY_Family=apply(amp[30:31,],2,sum),
  HSFX_Family=apply(amp[32:33,],2,sum),
  VCX_Family=apply(amp[34:37,],2,sum),
  DAZL=amp[38,],
  CDYL_Family=apply(amp[39:40,],2,sum),
  TSPYL2=amp[41,],
  RBMX_Family=apply(amp[42:44,],2,sum)
)
colnames(Amplicon_summary)==colnames(d)
d_amp<-rbind(d,Amplicon_summary[-c(12,14),])


family_Exp<-as.data.frame(t(Amplicon_summary),stringsAsFactors=FALSE,colClasses = "double")
#CDY~CDYL
#plot(Amplicon_summary[2,],Amplicon_summary[13,])
r=cor.test(family_Exp$CDY_Family,family_Exp$CDYL_Family,method = "spearman")
lrm<-lm(CDY_Family~CDYL_Family, family_Exp)
lrm.cc<-lrm$coefficients
CDY_ggO <- ggplot(family_Exp, aes(y=CDY_Family, x=CDYL_Family)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp, mapping=aes(y = CDY_Family, x = CDYL_Family), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("CDY=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "CDYL", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() + 
		theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5))+
		labs(title= "CDY")

#DAZ~DAZL
lrm<-lm(DAZ_Family~DAZL, family_Exp)
lrm.cc<-lrm$coefficients
r=cor.test(family_Exp$DAZ_Family,family_Exp$DAZL,method = "spearman")
DAZ_ggO <- ggplot(family_Exp, aes(y=DAZ_Family, x=DAZL)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp, mapping=aes(y = DAZ_Family, x = DAZL), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("DAZ=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "DAZL", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() + theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5))+labs(title= "DAZ")

#plot(Amplicon_summary[3,],Amplicon_summary[12,])


#HSFY~HSFX
HSFY_outlier<-which(family_Exp$HSFX_Family>1500)
family_Exp_HSFY<-family_Exp[-HSFY_outlier,]
lrm<-lm(HSFY_Family~HSFX_Family, family_Exp_HSFY)
lrm.cc<-lrm$coefficients
r=cor.test(family_Exp_HSFY$HSFY_Family,family_Exp_HSFY$HSFX_Family,method = "spearman")
HSFY_ggO <- ggplot(family_Exp_HSFY, aes(y=HSFY_Family, x=HSFX_Family)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp_HSFY, mapping=aes(y = HSFY_Family, x = HSFX_Family), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("HSFY=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "HSFX", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() + theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5))+labs(title= "HSFY")
#plot(Amplicon_summary[4,],Amplicon_summary[10,])
#RBMY~RBMX
r=cor.test(family_Exp$RBMY_Family,family_Exp$RBMX_Family,method = "spearman")
lrm<-lm(RBMY_Family~RBMX_Family, family_Exp)
lrm.cc<-lrm$coefficients
RBMY_ggO <- ggplot(family_Exp, aes(y=RBMY_Family, x=RBMX_Family)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp, mapping=aes(y = RBMY_Family, x = RBMX_Family), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("RBMY=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "RBMX", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() + theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5))+labs(title= "RBMY")
#plot(Amplicon_summary[6,],Amplicon_summary[15,])
#TSPY~TSPYL
lrm<-lm(TSPY_Family~TSPYL2, family_Exp)
lrm.cc<-lrm$coefficients
r=cor.test(family_Exp$TSPY_Family,family_Exp$TSPYL2,method = "spearman")
TSPY_ggO <- ggplot(family_Exp, aes(y=TSPY_Family, x=TSPYL2)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp, mapping=aes(y = TSPY_Family, x = TSPYL2), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("TSPY=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "TSPYL2", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() + theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5))+labs(title= "TSPY")
#plot(Amplicon_summary[7,],Amplicon_summary[14,])
#VCY~VCX
VCY_outlier<-which(family_Exp$VCY_Family==max(family_Exp$VCY_Family))
family_Exp_VCY<-family_Exp[-VCY_outlier,]
lrm<-lm(VCY_Family~VCX_Family,family_Exp_VCY )
lrm.cc<-lrm$coefficients
r=cor.test(family_Exp_VCY$VCY_Family,family_Exp_VCY$VCX_Family,method = "spearman")
VCY_ggO <- ggplot(family_Exp_VCY, aes(y=VCY_Family, x=VCX_Family)) + 
		geom_point(col=colorPalette[8]) +
		geom_smooth(data=family_Exp_VCY, mapping=aes(y = VCY_Family, x = VCX_Family), col="black",method=lm,se=FALSE) +
		#xlim(3,10)+ylim(5,11)+
		#annotate("text",x=-Inf,y=-Inf,size=4,vjust=-1, hjust=-0.1,label=paste("r=",round(r$estimate,digits = 2),"(P=",round(r$p.value,digits = 4),")",by=""))+
		annotate("text",x=-Inf,y=Inf,size=2.5,vjust=1, hjust=0,label=(eqn <- paste("VCY=", paste(round(lrm.cc[1],2), paste(round(lrm.cc[-1],2), "VCX", sep=" * ", collapse=" + "), sep=" + "))))+
		annotate("text",x=Inf,y=-Inf,size=2.5,vjust=0, hjust=1,label=paste("R2=",round(summary(lrm)$r.squared,digits = 2),by=""))+ 
		theme_bw() +labs(title= "VCY")
#plot(Amplicon_summary[8,],Amplicon_summary[11,])

pdf("Fig5B_GeneExpression_DOSAGE_protoXY_plot_newsize12.pdf")
grid.arrange(
CDY_ggO+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
DAZ_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
HSFY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
RBMY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
TSPY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
VCY_ggO+theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),ncol=3)
dev.off()


# CDY_LIST<-c("CDY_Family","CDYL","CDYL2")
# CDY_LIST_CHR<-c("Y","A","A")
# DAZ_LIST<-c("DAZ_Family","DAZL","BOLL")
# DAZ_LIST_CHR<-c("Y","A","A")
# HSFY_LIST_xfam<-c("HSFY_Family","HSFX_GF","HSFX1","HSFX2","HSFY1P1")
# HSFY_LIST_CHR_xfam<-c("Y","X","X","X","A")

# HSFY_LIST<-c("HSFY_Family","HSFX1","HSFX2","HSFY1P1")
# HSFY_LIST_CHR<-c("Y","X","X","A")
# RBMY_LIST<-c("RBMY_Family","RBMX","RBMX2","RBMXL3","RBMXL1","RBMXL2")
# RBMY_LIST_CHR<-c("Y","X","X","X","A","A")
# TSPY_LIST=c("TSPY_Family","TSPYL2","TSPYL1","TSPYL4","TSPYL5","TSPYL6")
# TSPY_LIST_CHR<-c("Y","X","A","A","A","A")
# VCY_LIST<-c("VCY_Family","VCX","VCX2","VCX3A","VCX3B")
# VCY_LIST_CHR<-c("Y","X","X","X","X")
# XKRY_LIST<-c("XKRY_Family","XKRX","XKR3")
# XKRY_LIST_CHR<-c("Y","X","A")


#Autosomal homolog
CDY_LIST<-c("CDY_Family","CDYL_Family") #"CDYL","CDYL2"
CDY_LIST_CHR<-c("Y","A")
DAZ_LIST<-c("DAZ_Family","DAZL")
DAZ_LIST_CHR<-c("Y","A")
#X homolog
HSFY_LIST<-c("HSFY_Family","HSFX_Family") #"HSFX1","HSFX2"
HSFY_LIST_CHR<-c("Y","X")
RBMY_LIST<-c("RBMY_Family","RBMX_Family") #"RBMX","RBMX2","RBMXL3"
RBMY_LIST_CHR<-c("Y","X")
TSPY_LIST=c("TSPY_Family","TSPYL2")
TSPY_LIST_CHR<-c("Y","X")
VCY_LIST<-c("VCY_Family","VCX_Family") #"VCX","VCX2","VCX3A","VCX3B"
VCY_LIST_CHR<-c("Y","X")

##XKRX is not expressed in testis. XKR3 is expressed in testis but not X homolog so not compared.
XKRY_LIST<-c("XKRY_Family","XKRX","XKR3")
XKRY_LIST_CHR<-c("Y","X","A")



generate_DF_GF<-function(GF_Homologs_List,GF_Homologs_List_CHR,GFname){
	GF_Homologs=NULL
	for (gene in GF_Homologs_List){
	   #print(gene)
	  #print(v[which(row.names(v)==gene),])
	  GF_Homologs=rbind(GF_Homologs,d_amp[which(row.names(d_amp)==gene),])
	}
	rownames(GF_Homologs)=GF_Homologs_List
	sam<-NULL
	data<- NULL
	gene<-NULL
	chr<-NULL
	for( i in 1:ncol(GF_Homologs)){
	  data<-c(data,GF_Homologs[,i])
	  sam<-c(sam,rep(colnames(GF_Homologs)[i],nrow(GF_Homologs)))
	  gene<-c(gene,rownames(GF_Homologs))
	  chr<-c(chr,GF_Homologs_List_CHR)
	}
	tab=cbind(data,sam,gene,chr)
	colnames(tab)=c("Expression","Sample","Gene","Chr")
	GF_df=as.data.frame(tab,stringsAsFactors = FALSE,row.names=F)
	GF_df$Expression=as.numeric(GF_df$Expression)
	cols <- c("A" = "#E69F00", "Y" = "#0078D7", "X" = "#009E73")  
	ggObject <- ggplot(GF_df, aes(y=Expression, x=Gene)) + 
		geom_boxplot(aes(color=Chr),outlier.shape=NA) +
		geom_point(aes(color=Chr),size=0.2) +
		scale_colour_manual(values = cols) +
		theme_bw() + 
		theme(                              
		axis.title.x = element_text(face="bold", color="black"),
		axis.title.y = element_text(face="bold", color="black"),
		plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),
		axis.text.x = element_text(angle = 45, hjust = 1,face = "bold"))+
		labs(title= paste(GFname)) +
		labs(x = "Genes") +
		labs(y = "Gene expression") +theme(legend.position="none")
	result <- list(data=GF_Homologs,plot=ggObject,ratio12=mean(GF_Homologs[1,])/mean(GF_Homologs[2,]),ratio21=mean(GF_Homologs[2,])/mean(GF_Homologs[1,]))
	return(result)
}

TSPY_LIST_summary<-generate_DF_GF(TSPY_LIST,TSPY_LIST_CHR,"TSPY")
DAZ_LIST_summary<-generate_DF_GF(DAZ_LIST,DAZ_LIST_CHR,"DAZ")
CDY_LIST_summary<-generate_DF_GF(CDY_LIST,CDY_LIST_CHR,"CDY")
HSFY_LIST_summary<-generate_DF_GF(HSFY_LIST,HSFY_LIST_CHR,"HSFY")
RBMY_LIST_summary<-generate_DF_GF(RBMY_LIST,RBMY_LIST_CHR,"RBMY")
VCY_LIST_summary<-generate_DF_GF(VCY_LIST,VCY_LIST_CHR,"VCY")
XKRY_LIST_summary<-generate_DF_GF(XKRY_LIST,XKRY_LIST_CHR,"XKRY")





pdf("Fig5A_GeneExpression_DOSAGE_protoXY_plot_newsize12.pdf")
grid.arrange(
CDY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),
DAZ_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),
HSFY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),
RBMY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),
TSPY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),
VCY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5),axis.text.x = element_text(angle = 20, hjust = 1)),ncol=2)
dev.off()

pdf("Fig5A_GeneExpression_DOSAGE_protoXY_plot_newsize8.pdf")
grid.arrange(
CDY_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
DAZ_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
HSFY_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
RBMY_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
TSPY_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
VCY_LIST_summary$plot+ theme_bw(base_size = 8)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),ncol=2)
dev.off()

pdf("Fig5_plot_newsize.pdf")
grid.arrange(
CDY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
DAZ_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
HSFY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
RBMY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
TSPY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
VCY_LIST_summary$plot+ theme_bw(base_size = 12)+theme(legend.position="none",aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
CDY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
DAZ_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
HSFY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
RBMY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
TSPY_ggO+ theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),
VCY_ggO+theme_bw(base_size = 12)+theme(aspect.ratio=1,plot.title = element_text(face="bold.italic", color = "black",, hjust = 0.5)),ncol=3)
dev.off()

ratioY<-cbind(
CDY_LIST_summary$ratio12,
DAZ_LIST_summary$ratio12,
HSFY_LIST_summary$ratio12,
RBMY_LIST_summary$ratio12,
TSPY_LIST_summary$ratio12,
VCY_LIST_summary$ratio12)

1-ratioY #Percentage difference between Y and homologs
#0.8928241 0.6295897 -0.4246934 0.6582184 0.7548914 0.7096196

