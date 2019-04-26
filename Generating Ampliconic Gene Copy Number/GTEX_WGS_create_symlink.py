#move files to data

import sys,os
import re

SRA_data="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/GTEx_WGS_SampleID_Run_Testis_Match.tab"
pathFolder="/nfs/secure/scratch6/ncbi/dbGaP-4719/sra"

SRAtoSampID=dict()
with open(SRA_data,"r") as i:
	for line in i:
		col=line.rstrip("\n").split("\t")
		#print col[0],col[1]
		SRAtoSampID[col[0]]=col[1]
		


moveto="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/"

Sample_list=os.listdir(pathFolder)
for id in Sample_list:
	p=os.path.join(pathFolder, id)
	if os.path.isfile(p):
		if id[:3]=="SRR":
			id_number=id.split("_")[0].split(".")[0]
			if id_number in SRAtoSampID:
				if id_number+"_Y.bam"==id:
					data_path=moveto+SRAtoSampID[id_number]+"/"
					if os.path.isdir(data_path):
						print "Already Created:"+SRAtoSampID[id_number]
						continue
					else:
						print "Checking BAM file"
						check_file=os.system("samtools quickcheck -v "+p+" && echo 'good' || echo 'fail'")
						print "Creating Folder:"
						print "mkdir "+data_path
						os.system("mkdir "+data_path)
						print "Copying BAM file:"
						print "ln -s "+p+" "+data_path
						os.system("ln -s "+p+" "+data_path)
						if os.path.isfile(p+".bai"):
							print "Copying Index file:"
							print "ln -s "+p+".bai "+data_path
							os.system("ln -s "+p+".bai "+data_path)
						else :
							print "ERROR: "+data_path+" Index file not found" 

