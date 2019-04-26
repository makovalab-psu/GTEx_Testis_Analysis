library(reshape2)
library(plyr)

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

##Obtaining family specific data in the form of list
geneFamilies <- split( Summary_AG , f = Summary_AG$Gene )



###MAIN TABLE

#Table 2 haplogroup information
##Haplogroup
table(YH_filtered$Haplogroup)
table(YH_filtered$SubHaplogroup)

#Table 3 ANOVA analysis
GeneLevel_HaplogroupvsGE_ANOVA<-function(genelevel){
	id=match(genelevel[,1],YH_filtered$Sample)
	#add haplogroup to gene level info
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	data_GF<-data_GF[(data_GF$Haplogroup %in% c("E","R","I","J")),]
	res.aov <- aov(Expression ~ Haplogroup, data = data_GF)
	out<-c(summary(res.aov),TukeyHSD(res.aov),pairwise.t.test(data_GF$Expression,data_GF$Haplogroup, p.adjust.method = "bonferroni"))
	return(out)
}

GeneLevel_HaplogroupvsCopyNumber_ANOVA<-function(genelevel){
	id=match(genelevel[,1],YH_filtered$Sample)
	#add haplogroup to gene level info
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	data_GF<-data_GF[(data_GF$Haplogroup %in% c("E","R","I","J")),]
	res.aov <- aov(CopyNumber ~ Haplogroup, data = data_GF)
	out<-c(summary(res.aov),TukeyHSD(res.aov),pairwise.t.test(data_GF$CopyNumber,data_GF$Haplogroup, p.adjust.method = "bonferroni"))
	return(out)
}

familylevelCopyNumbervsHaplo_test=llply(geneFamilies,GeneLevel_HaplogroupvsCopyNumber_ANOVA)
CNpvalueANOVA<-as.data.frame(llply(familylevelCopyNumbervsHaplo_test, function(x){c(x[[1]][["F value"]][[1]],x[[1]][["Pr(>F)"]][[1]])}))

familylevelGeneExpressionvsHaplo_test=llply(geneFamilies,GeneLevel_HaplogroupvsGE_ANOVA)
GEpvalueANOVA<-as.data.frame(llply(familylevelGeneExpressionvsHaplo_test, function(x){c(x[[1]][["F value"]][[1]],x[[1]][["Pr(>F)"]][[1]])}))

Table3<-t(rbind(CNpvalueANOVA,GEpvalueANOVA))
colnames(Table3)=c("CN.F-value","CN.P-value","GE.F-value","GE.P-value")
write.table(Table3,file="GTEx_ANOVA_Haplogroup.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


##############Supplementary Table

#Table S6. Median, standard deviation (SD) and range of copy number (CN, N=167) and gene expression (GE) values per ampliconic gene family (N=149).
medianAmpliconGE=apply(AGE_filtered,1,median)
sdAmpliconGE=apply(AGE_filtered,1,sd)
rangeAmpliconGE=apply(AGE_filtered,1,function(x){a<-range(x);paste(round(a[1]),round(a[2]),sep="-")})

medianAmpliconCN=apply(AmpliconCN_filtered,1,median)
sdAmpliconCN=apply(AmpliconCN_filtered,1,sd)
rangeAmpliconCN=apply(AmpliconCN_filtered,1,function(x){a<-range(x);paste(round(a[1],digit=2),round(a[2],digit=2),sep="-")})

TableS6=cbind(round(medianAmpliconCN,digit=2),round(sdAmpliconCN,digit=2),rangeAmpliconCN,round(medianAmpliconGE,digit=2),round(sdAmpliconGE,digit=2),rangeAmpliconGE)
colnames(TableS6)=c("CN.Median(N=167)","CN.SD","CN.Range","GE.Median(N=149)","GE.SD","GE.Range")
write.table(TableS6,file="GTEx_Median_SD_RANGE.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


#Table S7. Copy number and gene expression correlation values.
GeneLevel_all_CNvsGE_correlation<-function(genelevel){
  r=cor.test(genelevel$CopyNumber,genelevel$Expression,method = "spearman")
  val=c(round(r$estimate,digits = 2),round(r$p.value,digits = 4))
  #print(val)
  return(val)
}

GeneLevel_Yhap_CNvsGE_correlation<-function(genelevel,haplogroup){
  id=match(genelevel[,1],YH_filtered$Sample)
  #add haplogroup to gene level info
  data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
  data_GH<-data_GF[which(data_GF$SubHaplogroup==haplogroup),]
  r=cor.test(data_GH$CopyNumber,data_GH$Expression,method = "spearman")
  val=c(round(r$estimate,digits = 2),round(r$p.value,digits = 4))
  #print(val)
  return(val)
}

genelevel_corrAll=llply(geneFamilies,GeneLevel_all_CNvsGE_correlation)
genelevel_corrR1b=llply(geneFamilies,function(x){GeneLevel_Yhap_CNvsGE_correlation(x,haplogroup="R1b")})
genelevel_corrE1b=llply(geneFamilies,function(x){GeneLevel_Yhap_CNvsGE_correlation(x,haplogroup="E1b")})
genelevel_corrI1a=llply(geneFamilies,function(x){GeneLevel_Yhap_CNvsGE_correlation(x,haplogroup="I1a")})

#Table S8. P-values from permutation tests for copy number differences between haplogroup pairs.
one.test <- function(x,y,h) {
xshuffeled<-sample(x)
mean(y[xshuffeled==h[1]])-mean(y[xshuffeled==h[2]])
}
GeneLevel_HaplogroupvsCopyNumberPermutationTest<-function(genelevel,haplogroup_list){
	set.seed(9)
	id=match(genelevel[,1],YH_filtered$Sample)
	#add haplogroup to gene level info
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	data_GF<-data_GF[(data_GF$Haplogroup %in% haplogroup_list),]
	true<-mean(data_GF$CopyNumber[data_GF$Haplogroup==haplogroup_list[1]])-mean(data_GF$CopyNumber[data_GF$Haplogroup==haplogroup_list[2]])
	permutation <- replicate(1000000, one.test(data_GF$Haplogroup,data_GF$CopyNumber ,haplogroup_list))
	pvalue<-mean(abs(true) < abs(permutation))
	return(pvalue)
	}



genelevel_permutationER=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("E","R"))})
genelevel_permutationEI=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("E","I"))})
genelevel_permutationEJ=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("E","J"))})
genelevel_permutationJR=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("J","R"))})
genelevel_permutationJI=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("J","I"))})
genelevel_permutationIR=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsCopyNumberPermutationTest(x,c("I","R"))})

#Table S9. P-values from permutation test for gene expression differences between haplogroup pairs. 

GeneLevel_HaplogroupvsGeneExpressionPermutationTest<-function(genelevel,haplogroup_list){
	set.seed(9)
	id=match(genelevel[,1],YH_filtered$Sample)
	#add haplogroup to gene level info
	data_GF<-cbind(genelevel,YH_filtered[id,c(2,3)])
	data_GF<-data_GF[(data_GF$Haplogroup %in% haplogroup_list),]
	true<-mean(data_GF$Expression[data_GF$Haplogroup==haplogroup_list[1]])-mean(data_GF$Expression[data_GF$Haplogroup==haplogroup_list[2]])
	permutation <- replicate(1000000, one.test(data_GF$Haplogroup,data_GF$Expression ,haplogroup_list))
	pvalue<-mean(abs(true) < abs(permutation))
	return(pvalue)
	}


genelevelExp_permutationER=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("E","R"))})
genelevelExp_permutationEI=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("E","I"))})
genelevelExp_permutationEJ=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("E","J"))})
genelevelExp_permutationJR=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("J","R"))})
genelevelExp_permutationJI=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("J","I"))})
genelevelExp_permutationIR=llply(geneFamilies,function(x){GeneLevel_HaplogroupvsGeneExpressionPermutationTest(x,c("I","R"))})


#Table S10. Correlation between gene expression and age.

Age_CNvsGE_correlation<-function(genelevel){
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
	val=c(round(r$estimate,digits = 2),round(r$p.value,digits = 4))
	#print(val)
	return(val)
}  

Age_CNvsGE_Haplogroup_correlation<-function(genelevel,haplogroup){
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
	val=c(round(r$estimate,digits = 2),round(r$p.value,digits = 4))
	#print(val)
	return(val)
}

genelevelAge_corrAll=llply(geneFamilies,Age_CNvsGE_correlation)
genelevelAge_corrR1b=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup_correlation(x,haplogroup="R1b")})
genelevelAge_corrE1b=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup_correlation(x,haplogroup="E1b")})
genelevelAge_corrI1a=llply(geneFamilies,function(x){Age_CNvsGE_Haplogroup_correlation(x,haplogroup="I1a")})
