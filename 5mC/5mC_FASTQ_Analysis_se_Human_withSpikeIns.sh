#!/bin/bash                                                                                                                                                  
                                       
#PBS -l nodes=1:ppn=6                                                                                                                                        
                                     
#PBS -l walltime=3:00:00:00


#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/Plate36_H9_scMAT_All_SpikeInTest/P36_H9_MouseBrainSpikeIn_RNA_5hmC_5mC_L1_All"




#Standard Usage, no input required
BARCODES="/home/achialastri/perlscripts/MspJI/mspj1"
#GENOME="/home/sdey/genomes/human_gene_models/hg19.fa"
GENOME="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage_pUC19_mm10/hg19_Zymo_LambdaPhage_pUC19_mm10.fa"
PERL_DIR="/home/achialastri/perlscripts/MspJI"
ABASI_BC="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"
CEL_BARCODES="/home/achialastri/perlscripts/mRNA_Mapping/cel-seq_barcodes.csv"


#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2

intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-MSPJI
OUT_NAME_R2=$intermediateR2-MSPJI

#Cat Fastq Files
#cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1
cp $START_DIR/*R1.fastq $START_DIR/$FASTQ_R1


#Deduplicate the fastq file
/home/achialastri/BisulfiteTools/bbmap/clumpify.sh in=$START_DIR/$FASTQ_R1 out=$START_DIR/$intermediateR1-dedup.fastq dedupe subs=0
rm $START_DIR/$FASTQ_R1
mv $START_DIR/$intermediateR1-dedup.fastq $START_DIR/$FASTQ_R1

#Pull out only MspJI Lines (Set --scMATseq == 1 if done with RNA, set --scTHseq == 1 if done with ABASI
perl $PERL_DIR/ExtractingMSPJIReads_UserInput_NoBarcodeCollisions.pl --FASTQ_R1 $START_DIR/$FASTQ_R1 --CELSEQ_BC $CEL_BARCODES --scMATseq 1 --MSPJI_BC $BARCODES --scTHseq 1 --ABASI_BC $ABASI_BC

#Mapping
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -q 0 -n 0.04 -k 2 -l 200 -t 6 -B 11 $GENOME $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME_R1.sai


/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n 100 $GENOME $START_DIR/$OUT_NAME_R1.sai $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME-se.sam

#Analysis prior to PCR can comment out if not interested in pre PCR information
perl $PERL_DIR/process_scmspjiwithQC_MouseAndHuman.pl $GENOME $START_DIR/$OUT_NAME-se.sam $BARCODES
/home/cwangsanuwat/src/samtools/samtools flagstat $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.flagstat


#Same as Bam File
/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.bam

#Separate out 5mC Types and Finialize location
perl $PERL_DIR/Make_Final_Version2.pl $START_DIR/$OUT_NAME-se.faba

#Counts in Faba
perl /home/achialastri/perlscripts/Reads_Per_Cell_Count.pl $START_DIR/$OUT_NAME-se.faba
perl /home/achialastri/perlscripts/Reads_Per_Cell_Count_txt.pl $START_DIR/$OUT_NAME-se_Full5mC_CG_Final_Rmdup.txt

#Remove excess files created
rm $START_DIR/$OUT_NAME-se.sam

rm $START_DIR/$OUT_NAME_R1.sai

rm $START_DIR/$OUT_NAME_R1.fastq 
rm $START_DIR/$OUT_NAME_R2.fastq 
 
rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2

#copy all important files to a single location for that run
FABA="_FABA"
FLAGSTAT="_FLAGSTAT"
RABA="_RABA"


mkdir $PAST_DIR/$RUN_NAME$FABA
cp $START_DIR/$OUT_NAME-se.faba $PAST_DIR/$RUN_NAME$FABA

mkdir $PAST_DIR/$RUN_NAME$FLAGSTAT
cp $START_DIR/$OUT_NAME-se.flagstat $PAST_DIR/$RUN_NAME$FLAGSTAT

mkdir $PAST_DIR/$RUN_NAME$RABA
cp $START_DIR/$OUT_NAME-se.raba $PAST_DIR/$RUN_NAME$RABA




