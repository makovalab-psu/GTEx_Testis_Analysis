#!/bin/bash
#
#SBATCH -o /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/Analysis/GTEX-ZDTT/SRR1466260.out
#SBATCH -D /nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/RNASeq/Analysis/GTEX-ZDTT
#SBATCH -J GTEX-ZDTT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=7000
#SBATCH --mail-type=end
#SBATCH --time=28:00:00

srun /nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/TESTIS/RNASeq/src/kallisto_linux-v0.43.0/kallisto quant --index=/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/TESTIS/RNASeq/ref/RNA/refMrna_kallisto_index --output-dir=kallisto_ZDTT --threads=4 --bootstrap-samples=100 --seed=9 --bias SRR1466260_1.fastq.gz SRR1466260_2.fastq.gz
