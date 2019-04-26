import os
import time;

#Create slurm jobs to process the files.


SRA_data="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/GTEx_WGS_SampleID_Run_Testis_Match.tab"

SampIDtoSRA=dict()
with open(SRA_data,"r") as i:
	for line in i:
		col=line.rstrip("\n").split("\t")
		SampIDtoSRA[col[1]]=col[0]




path_data="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/" #path should end with "/"
Sample_list=os.listdir(path_data)

work_dir="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/analysis/ProcessBAM"

status_file="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/analysis/ProcessBAM/Status_Folder_Creation.txt"
sta=open(status_file,"a")
localtime = time.localtime(time.time())
print "Local current time :", localtime
sta.write("\n\nLocal current time :"+str(localtime)+"\n")
for sid in Sample_list:
	p=os.path.join(path_data, sid)
	if os.path.isdir(p):
		dest_folder=work_dir+"/"+sid
		if os.path.isdir(dest_folder):
		#	print sid+"\tFolder already exist"
		#else:
			print sid+"\tCreated"
			sta.write(sid+"\tCreated\n")
			#os.system("mkdir "+dest_folder)
			BAM_PATH=p+"/"+SampIDtoSRA[sid]+"_Y.bam"
			#print BAM_PATH
			#Generate SLURM job			
			#sl=open(dest_folder+"/slurm_Amplicon_Haplo.sh","w")
			sl=open(dest_folder+"/slurm_AmpliCoNE.sh","w")
			sl.write("#!/bin/bash\n")
			sl.write("#\n")
			sl.write("#SBATCH -o "+work_dir+"/"+sid+"/"+sid+"_slurm_AmpliCoNE.out\n")
			sl.write("#SBATCH -D "+work_dir+"/"+sid+"\n")
			sl.write("#SBATCH -J "+sid+"\n")
			sl.write("#SBATCH --nodes=1\n")
			sl.write("#SBATCH --ntasks=1\n")
			sl.write("#SBATCH --cpus-per-task=10\n")
			sl.write("#SBATCH --mem=83000\n")
			sl.write("#SBATCH --mail-type=end\n")
			#sl.write("#SBATCH --mail-user=<youreamilid>\n")
			sl.write("#SBATCH --time=24:00:00\n\n")
			sl.write("\ntime srun python /galaxy/home/rxv923/scripts/cn/AmpliCoNE.py --BAM "+SampIDtoSRA[sid]+"_Ychr_sorted_RMDup.bam --CHR Y\n")
			#sl.write("\n#Steps too extract reads\nsrun samtools view -hb -f 3 -F 12 "+BAM_PATH+" > Parsed_Ychr_pairedReads.bam\nsrun samtools sort -n Parsed_Ychr_pairedReads.bam > Parsed_Ychr_pairedReads_sotredbyPOS.bam\nsrun bamToFastq -i Parsed_Ychr_pairedReads_sotredbyPOS.bam -fq "+SampIDtoSRA[sid]+"_Ychr_1.fastq -fq2 "+SampIDtoSRA[sid]+"_Ychr_2.fastq\n\n#Amplicone Steps\nsrun sh /galaxy/home/rxv923/scripts/cn/GTEx_pipeline.sh BWA "+SampIDtoSRA[sid]+"_Ychr "+sid+"\ntime srun python /galaxy/home/rxv923/scripts/cn/AmpliCoNE.py --BAM "+SampIDtoSRA[sid]+"_Ychr_sorted_RMDup.bam --CHR Y\n")
			#sl.write("\n#Haplogroup Steps\nsrun samtools mpileup -f /nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/Homo_sapiens_assembly19.fasta "+BAM_PATH+" -r \"Y:2655100-28771000\" >"+SampIDtoSRA[sid]+"_hg19_chrY_2655100_28771000.pileup\nsrun python /galaxy/home/rxv923/scripts/Haplo/make_input.py "+SampIDtoSRA[sid]+"_hg19_chrY_2655100_28771000.pileup "+SampIDtoSRA[sid]+"_hg19.genos.txt\nsrun python /galaxy/home/rxv923/src/yhaplo/callHaplogroups.py -i "+SampIDtoSRA[sid]+"_hg19.genos.txt")
	else:
		print sid+" is not a Directory"

sta.close()
