import sys,os
import re

###########################################
#Usage : python GTEX_WGS_parse_outputfiles.py
#Parse the output of AmpliCoNE and Yhaplo for 170 files.
#
#"pathFolder" is folder where the outputs of AmpliCoNE and Yhaplo are stored for each sample.
#	structure: /pathFolder/
#					<SAMPLE 1>/
#						<SAMPLE ID>__sorted_RMDup.bamAmpliconic_Summary.txt --> AmpliCoNE output Gene Family summary
#						<SAMPLE ID>_sorted_RMDup.bamXDG_CopyNumber.txt		--> AmpliCoNE output XDG summary
#						/output/haplogroups.<SAMPLE ID>.txt					--> Yhaplo output   
#						.*Other intermediate file : FASTQ, BAM, BAI, pileup, etc. 
#						
#					<SAMPLE 2>/
#						<SAMPLE ID>__sorted_RMDup.bamAmpliconic_Summary.txt --> AmpliCoNE output Gene Family summary
#						<SAMPLE ID>_sorted_RMDup.bamXDG_CopyNumber.txt		--> AmpliCoNE output XDG summary
#						/output/haplogroups.<SAMPLE ID>.txt					--> Yhaplo output   
#						.*Other intermediate file : FASTQ, BAM, BAI, pileup, etc.
							
#Final output: Folder "moveto" with all the output files and file "Haplogroup_summary" with haplogroups info for all the samples.

###########################################

#Path to where the BAM files are located
pathFolder="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/analysis/ProcessBAM"

#File with SRA ID and subject id
#Run     submitted_subject_id
#SRR2166358      GTEX-PLZ5
#SRR2166862      GTEX-PLZ6
#SRR2165091      GTEX-TKQ2
SRA_data="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/GTEx_WGS_SampleID_Run_Testis_Match.tab"

SampIDtoSRA=dict()
with open(SRA_data,"r") as i:
	for line in i:
		col=line.rstrip("\n").split("\t")
		SampIDtoSRA[col[1]]=col[0]

#Parse AmpliCoNE output files and move them to a single folder.
moveto="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output/"
#File into which the haplogroup summary will be printed
Haplogroup_summary="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/results/GTEx_WGS_ParseYchr/AmpliCoNE_output/GTEX_Yhaplogroups_YchrBAM.tab"
hg = open(Haplogroup_summary, 'w')
Sample_list=os.listdir(pathFolder)
for id in Sample_list:
	p=os.path.join(pathFolder, id)
	if os.path.isdir(p):
		#print moveto+id
		
		for path, subdirs, files in os.walk(p):
			for name in files:
				if bool(re.search("Ampliconic_Summary.txt",name)): #Name of the output file can change in the future, update it here.
					#print os.path.join(path, name)
					#print "cp "+os.path.join(path, name)+" "+moveto+id
					print name
					os.system("mkdir "+moveto+id) #create the folder if it does not exist
					os.system("cp "+os.path.join(path, name)+" "+moveto+id) #Move the file to single location.
				elif bool(re.search("XDG_CopyNumber.txt",name)):  #Name of the output file can change in the future, update it here.
					#print os.path.join(path, name)
					print name
					#os.system("mkdir "+moveto+id)
					os.system("cp "+os.path.join(path, name)+" "+moveto+id)
				elif bool(re.search("haplogroups.",name)): #Name of the output file can change in the future, update it here.
					#print os.path.join(path, name)
					print name
					#os.system("mkdir "+moveto+id)
					#os.system("cp "+os.path.join(path, name)+" "+moveto+id)
					f = open(os.path.join(path, name), 'r')
					output_haplo=f.read().split()
					hg.write(id+"\t"+SampIDtoSRA[id]+"\t"+output_haplo[1]+"\t"+output_haplo[2]+"\t"+output_haplo[3]+"\n")


hg.close()
