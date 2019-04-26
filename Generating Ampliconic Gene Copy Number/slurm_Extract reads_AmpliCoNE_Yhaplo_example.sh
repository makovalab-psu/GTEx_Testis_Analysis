#!/bin/bash
#
#SBATCH -o /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/analysis/ProcessBAM/GTEX-ZDTT/GTEX-ZDTT_slurm.out
#SBATCH -D /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/analysis/ProcessBAM/GTEX-ZDTT
#SBATCH -J GTEX-ZDTT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=83000
#SBATCH --mail-type=end
#SBATCH --time=24:00:00


#Steps too extract reads
srun samtools view -hb -f 3 -F 12 /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/GTEX-ZDTT/SRR3481969_Y.bam > Parsed_Ychr_pairedReads.bam
srun samtools sort -n Parsed_Ychr_pairedReads.bam > Parsed_Ychr_pairedReads_sotredbyPOS.bam
srun bamToFastq -i Parsed_Ychr_pairedReads_sotredbyPOS.bam -fq SRR3481969_Ychr_1.fastq -fq2 SRR3481969_Ychr_2.fastq

#Amplicone Steps
srun sh /galaxy/home/rxv923/scripts/cn/GTEx_pipeline.sh BWA SRR3481969_Ychr GTEX-ZDTT
time srun python /galaxy/home/rxv923/scripts/cn/AmpliCoNE.py --BAM SRR3481969_Ychr_sorted_RMDup.bam --CHR Y

#Haplogroup Steps
srun samtools mpileup -f /nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/Homo_sapiens_assembly19.fasta /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/data/BAM_GTEX_Ychr/GTEX-ZDTT/SRR3481969_Y.bam -r "Y:2655100-28771000" >SRR3481969_hg19_chrY_2655100_28771000.pileup
srun python /galaxy/home/rxv923/scripts/Haplo/make_input.py SRR3481969_hg19_chrY_2655100_28771000.pileup SRR3481969_hg19.genos.txt
srun python /galaxy/home/rxv923/src/yhaplo/callHaplogroups.py -i SRR3481969_hg19.genos.txt
