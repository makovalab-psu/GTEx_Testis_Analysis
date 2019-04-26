#!usr/bin/sh

# bail out on errors
set -ue

STEP=0
LEVEL=$1

SAMPLE=$2	#SRR2165027
ID=$3			#GTEX-N7MS
THREADS=10
#the path to the reference file, the folder should have the BWA index and FAI file in it along with the reference
REF=/galaxy/home/rxv923/refs/Hg/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa

#Set the exact path to where the binaries are for the tools 
BWA=/galaxy/home/rxv923/src/bwa-0.7.15/bwa
PICARD=/galaxy/home/rxv923/src/picard-master/dist/picard.jar
SAMTOOLS=/galaxy/home/rxv923/src/samtools-1.3.1/bin/samtools

#Folder where there is enough space to store intermediate files.
TEMPDIR=/nfs/secure/scratch6/nekrut_gtex/rxv923/GTEX/TESTIS/WGS/TEST_dbGAP/



case $LEVEL in
        "FASTQC")
                STEP=1;;
        "FIXFASTQ")
                STEP=2;;
        "BWA")
                STEP=3;;
        "SORTBAM")
                STEP=4;;
        "INDEXBAM")
                STEP=5;;
        "RMDUP")
                STEP=6;;
        "INDEXBAM2")
                STEP=7;;
        *);;
esac


#fastq-dump -I --split-files SRR2167248.sra
if [ $STEP -le 1 ]
then
        echo "Generating Quality report for the FASTQ files"
        fastqc ${SAMPLE}_1.fastq
        fastqc ${SAMPLE}_2.fastq
fi
# if [ $STEP -le 2 ]
# then
        # echo "Fixing the FastQ file by renaming the last two letters in the read name from .1/2 to <space>1/2"
        # python /galaxy/home/rxv923/scripts/flyscripts/fix_fastq.py ${SAMPLE}_1.fastq
        # python /galaxy/home/rxv923/scripts/flyscripts/fix_fastq.py ${SAMPLE}_2.fastq
# fi
if [ $STEP -le 3 ]
then
        echo "Align the reads using BWA MEM"
        RG="@RG\tID:${SAMPLE}\tSM:${ID}"
        echo $RG
        $BWA mem -t $THREADS -R $RG $REF ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq > ${SAMPLE}.sam
fi
if [ $STEP -le 4 ]
then
        echo "Sorting BAM file by location"
        java -Xmx4g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=${SAMPLE}.sam OUTPUT=${SAMPLE}_sorted.bam TMP_DIR=$TEMPDIR
fi
if [ $STEP -le 5 ]
then
        echo "Indexing BAM file"
        $SAMTOOLS index ${SAMPLE}_sorted.bam
fi
if [ $STEP -le 6 ]
then
        echo "Removing Duplicate Reads"
        java  -Xmx10g -XX:-UseGCOverheadLimit -jar  $PICARD MarkDuplicates REMOVE_DUPLICATES=TRUE INPUT=${SAMPLE}_sorted.bam OUTPUT=${SAMPLE}_sorted_RMDup.bam METRICS_FILE=${SAMPLE}_duplicate_metrics.txt TMP_DIR=$TEMPDIR  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
fi
if [ $STEP -le 7 ]
then
        echo "Indexing BAM file"
        $SAMTOOLS index ${SAMPLE}_sorted_RMDup.bam
fi
