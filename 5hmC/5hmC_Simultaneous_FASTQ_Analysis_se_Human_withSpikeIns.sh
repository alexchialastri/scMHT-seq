#!/bin/bash                                                                                                                                                  
                                       
#PBS -l nodes=1:ppn=6                                                                                                                                        
                                     
#PBS -l walltime=3:00:00:00


#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/Plate36_H9_scMAT_All_SpikeInTest/P36_H9_MouseBrainSpikeIn_RNA_5hmC_5mC_L1_All"




#Standard Usage, no input required
BARCODES="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"
#GENOME="/home/sdey/genomes/human_gene_models/hg19.fa"
GENOME="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage_pUC19_mm10/hg19_Zymo_LambdaPhage_pUC19_mm10.fa"
PERL_DIR="/home/achialastri/perlscripts/5hmC"
MSPJIBARCODES="/home/achialastri/perlscripts/MspJI/mspj1"
CEL_BARCODES="/home/achialastri/perlscripts/mRNA_Mapping/cel-seq_barcodes.csv"


#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}-ABA
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2

intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-ABA
OUT_NAME_R2=$intermediateR2-ABA

#Cat Fastq Files
#cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1
cp $START_DIR/*R1.fastq $START_DIR/$FASTQ_R1


#Deduplicate the fastq file
/home/achialastri/BisulfiteTools/bbmap/clumpify.sh in=$START_DIR/$FASTQ_R1 out=$START_DIR/$intermediateR1-dedup.fastq dedupe subs=0
rm $START_DIR/$FASTQ_R1
mv $START_DIR/$intermediateR1-dedup.fastq $START_DIR/$FASTQ_R1



#Order Matters for Arguments
#perl $PERL_DIR/ExtractingAbaReads_96BC_SimultanousWith5mC_UserInput.pl $START_DIR $OUT_NAME_R1 $OUT_NAME_R2 $FASTQ_R1 $FASTQ_R2 $BARCODES $MSPJIBARCODES
#Pull out only AbaSI Lines (Set --scMHseq == 1 if done with mspji, set --scTHseq == 1 if done with RNA
perl $PERL_DIR/ExtractingABASIReads_UserInput_NoBarcodeCollisions.pl --FASTQ_R1 $START_DIR/$FASTQ_R1 --CELSEQ_BC $CEL_BARCODES --scMHseq 1 --MSPJI_BC $MSPJIBARCODES --scTHseq 1 --ABASI_BC $BARCODES

#Mapping
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -q 0 -n 0.04 -k 2 -l 200 -t 6 -B 6 $GENOME $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME_R1.sai


/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n 100 $GENOME $START_DIR/$OUT_NAME_R1.sai $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME-se.sam

#call 5hmC sites
perl $PERL_DIR/process_scaba_AnyInput_MouseAndHuman.pl $GENOME $START_DIR/$OUT_NAME-se.sam  $PERL_DIR/aba_barcodes.csv

#mapping info
/home/cwangsanuwat/src/samtools/samtools flagstat $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se-ABA.flagstat


#Same as Bam File
/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se-ABA.bam



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
cp $START_DIR/$OUT_NAME-se-ABA.faba $PAST_DIR/$RUN_NAME$FABA

mkdir $PAST_DIR/$RUN_NAME$FLAGSTAT
cp $START_DIR/$OUT_NAME-se-ABA.flagstat $PAST_DIR/$RUN_NAME$FLAGSTAT

mkdir $PAST_DIR/$RUN_NAME$RABA
cp $START_DIR/$OUT_NAME-se-ABA.raba $PAST_DIR/$RUN_NAME$RABA




