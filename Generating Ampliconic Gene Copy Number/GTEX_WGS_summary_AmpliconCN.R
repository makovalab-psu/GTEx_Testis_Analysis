library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggsignif)

#Create a summary file for AmpliCoNE output, Table where row is gene family and column is sample. Two files will be generated one has Ampliconic gene family summary and other 
#Usage: Rscript GTEX_WGS_summary_AmpliconCN.R 
#samples_path is the output folder (moveto) from GTEX_WGS_parse_outputfiles.py
#Haplogroup is the output file (Haplogroup_summary) from GTEX_WGS_parse_outputfiles.py

samples_path<-list.dirs(path = "/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output", recursive=FALSE) #no / at the end of the path
#This is same as samples_path. Used to get list of files inside the folder, for which the output is parsed and summarized.
samples<-list.dirs(path = "/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output", recursive=FALSE,full.names = FALSE)
Haplogroup<-"/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output/GTEX_Yhaplogroups_YchrBAM.tab"

Yhapl<-read.table(Haplogroup, sep="\t", stringsAsFactors=FALSE,header=FALSE)
Yhapl<-cbind(Yhapl,substr(Yhapl[,5],1,1))
colnames(Yhapl)[c(5,6)]<-c("YHaplo_True","Haplogroups")
AmpCN<-NULL
XdgCN<-NULL

for( i in seq(length(samples_path))){
	print(samples_path[i])
	print(samples[i])
	output_files<-list.files(samples_path[i],all.files = FALSE,full.names = FALSE)
	temp_fn_a<-paste(samples_path[i],output_files[grep("Ampliconic_Summary.txt",output_files)],sep="/") #The name of the AampliCoNE output might change, update here.
	temp_acn<-read.table(temp_fn_a, sep="\t", stringsAsFactors=FALSE,header=TRUE)
	AmpCN<-cbind(AmpCN,temp_acn[,2])
	colnames(AmpCN)[i]=samples[i]
	
	temp_fn_x<-paste(samples_path[i],output_files[grep("XDG_CopyNumber.txt",output_files)],sep="/")  #The name of the AampliCoNE output might change, update here.
	temp_xcn<-read.table(temp_fn_x, sep="\t", stringsAsFactors=FALSE,header=TRUE)
	XdgCN<-cbind(XdgCN,temp_xcn[,2])
	colnames(XdgCN)[i]=samples[i]
	
}
AmpCN_file=cbind(temp_acn[,1],AmpCN)
colnames(AmpCN_file)[1]="GeneFamily"
rownames(AmpCN)=temp_acn[,1]
XdgCN_file=cbind(temp_xcn[,1],XdgCN)
colnames(XdgCN_file)[1]="GeneFamily"
rownames(XdgCN)=temp_xcn[,1]

write.table(AmpCN_file,file="GTEx_CopyNumber_Amplicon.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(XdgCN_file,file="GTEx_CopyNumber_XDG.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


