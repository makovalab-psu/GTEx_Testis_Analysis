import os


#This code reads a file, parses the sample id and SRR id, check if file exists in the given folder and creates a symlink and folder and move the files there.
#path_data="/nfs/secure/scratch4/boris/amelia/36828/reads/fqgz/" #path should end with "/" ##Folder where Boris downloaded the data
path_data="/nfs/secure/scratch6/ncbi/dbGaP-4719/sra/" #path should end with "/" ##Folder where Rico downloaded the data
work_dir="/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/Analysis/" #path should end with "/" 

#with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR.tab",'r') as input: #####Samples analyzed in second iteration
#with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR_First23.tab",'r') as input: #####Samples analyzed first iteration
with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR_recoDownloadSet.tab",'r') as input: #####Samples analyzed in second iteration
	for line in input:
		 #os.system('gzip '+line)
		col=line.rstrip('\n').split("\t")
		file=os.path.isfile(path_data+col[0]+"_1.fastq.gz")
		if file:
			if os.path.isdir(work_dir+col[1]):
				print "Folder Exists: ",col
			else:
				print "Created Link: ",col[1]
				os.system("ln -s "+path_data+col[0]+"_1.fastq.gz .")
				os.system("ln -s "+path_data+col[0]+"_2.fastq.gz .")
				os.system("mkdir "+work_dir+col[1])
				os.system("mv "+col[0]+"_* "+work_dir+col[1])
				sl=open(work_dir+col[1]+"/slurm_kallisto.sh","w")
				sl.write("#!/bin/bash\n")
				sl.write("#\n")
				sl.write("#SBATCH -o "+work_dir+col[1]+"/"+col[0]+".out\n")
				sl.write("#SBATCH -D "+work_dir+col[1]+"\n")
				sl.write("#SBATCH -J "+col[1]+"\n")
				sl.write("#SBATCH --nodes=1\n")
				sl.write("#SBATCH --ntasks=1\n")
				sl.write("#SBATCH --cpus-per-task=5\n")
				sl.write("#SBATCH --mem=7000\n")
				sl.write("#SBATCH --mail-type=end\n")
				#sl.write("#SBATCH --mail-user=<youreamilid>\n")
				sl.write("#SBATCH --time=28:00:00\n\n")
				sl.write("srun /nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/TESTIS/RNASeq/src/kallisto_linux-v0.43.0/kallisto quant --index=/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/TESTIS/RNASeq/ref/RNA/refMrna_kallisto_index --output-dir=kallisto_"+col[1][5:]+" --threads=4 --bootstrap-samples=100 --seed=9 --bias "+col[0]+"_1.fastq.gz "+col[0]+"_2.fastq.gz\n")
		else:
				print "ERROR: FASTQ file is not found", col
