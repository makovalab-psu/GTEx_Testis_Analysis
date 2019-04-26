import os


#This code reads a file, parses the sample id and SRR id, check if file exists in the given folder and creates a symlink and folder and move the files there.
#path_data="/nfs/secure/scratch4/boris/amelia/36828/reads/fqgz/" #path should end with "/" ##Folder where Boris downloaded the data
path_data="/nfs/secure/scratch6/ncbi/dbGaP-4719/sra/" #path should end with "/" ##Folder where Rico downloaded the data


#with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR.tab",'r') as input: #####Samples analyzed in second iteration
#with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR_First23.tab",'r') as input: #####Samples analyzed first iteration
with open("/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/sampleID_SRR_recoDownloadSet.tab",'r') as input: #####Samples analyzed in second iteration
        for line in input:
                #os.system('gzip '+line)
                col=line.rstrip('\n').split("\t")
                file=os.path.isfile(path_data+col[0]+"_1.fastq.gz")
                if file:
                        print "Created Link: ",col
                        os.system("ln -s "+path_data+col[0]+"_1.fastq.gz .")
                        os.system("ln -s "+path_data+col[0]+"_2.fastq.gz .")
                        os.system("mkdir /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/"+col[1])
                        os.system("mv "+col[0]+"_* /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/data/"+col[1])
                else:
                        print "ERROR: File not found", col
