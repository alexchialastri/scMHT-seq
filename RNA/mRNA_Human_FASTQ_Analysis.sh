#!/bin/bash                                                                                                                                                                                         
#PBS -l nodes=1:ppn=6                                                                                                                                                                             
#PBS -l walltime=24:00:00


#User Input
START_DIR="/home/achialastri/PGCLC_Extended_UCLA1_WT/PGCLCF6_L2R_UCLA1_D4"



#Standard Usage per person
CEL_BARCODES="/home/achialastri/perlscripts/mRNA_Mapping/cel-seq_barcodes.csv"
TRANSCRIPTOME="/home/achialastri/perlscripts/mRNA_Mapping/hg19_RefSeq_clean_ERCC92.fa"
IsoToGene="/home/achialastri/perlscripts/mRNA_Mapping/hg19_RefSeq2gene.tsv"
PERL_DIR="/home/achialastri/perlscripts/mRNA_Mapping"
HOME_DIR="/home/achialastri"
ABASI_BC="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"
MSPJI_BC="/home/achialastri/perlscripts/MspJI/mspj1"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2
LOG=${FASTQ_R1%??????}
LOG2=${FASTQ_R2%??????}


UMILength=6

#Cat Fastq Files
#cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1.all
#cat $START_DIR/*L001_R2* $START_DIR/*L002_R2* $START_DIR/*L003_R2* $START_DIR/*L004_R2* > $START_DIR/$FASTQ_R2.all
cp $START_DIR/*_R1.fastq $START_DIR/$FASTQ_R1.all
cp $START_DIR/*_R2.fastq $START_DIR/$FASTQ_R2.all

#Pull out only Celseq Lines (Set --scMATseq == 1 if done with MSPJI, set --scTHseq == 1 if done with ABASI
perl $PERL_DIR/Extracting_mRNA_Reads_96BC_UserInput_NoBarcodeCollisions.pl --StartingDir $START_DIR --OUTNAME_R1 $LOG --OUTNAME_R2 $LOG2 --FASTQ_R1 $FASTQ_R1.all --FASTQ_R2 $FASTQ_R2.all --CELSEQ_BC $CEL_BARCODES --scMATseq 1 --MSPJI_BC $MSPJI_BC --scTHseq 1 --ABASI_BC $ABASI_BC

perl $PERL_DIR/do_mappings_strand_Deylab.pl -r=$TRANSCRIPTOME -f1=$START_DIR/$FASTQ_R1 -f2=$START_DIR/$FASTQ_R2 -out=$OUT_NAME -outdir=$START_DIR/$OUT_NAME -t=6 -cel=1 -bar=$CEL_BARCODES -fstr=1 -rb=1 -rb_len=$UMILength -nsam=0 > $LOG.log1 2> $LOG.log2

mv $HOME_DIR/$LOG.log1 $START_DIR/$OUT_NAME
mv $HOME_DIR/$LOG.log2 $START_DIR/$OUT_NAME

perl $PERL_DIR/iso2gene_rb.pl -m=$IsoToGene -in=$START_DIR/$OUT_NAME/$OUT_NAME.cout.csv -out=$OUT_NAME.gene.cout.csv

mv $HOME_DIR/$OUT_NAME.gene.cout.csv $START_DIR

perl $PERL_DIR/extract_counts_rb.pl -bl=$UMILength -in=$START_DIR/$OUT_NAME.gene.cout.csv -outc=$OUT_NAME-READ_COUNTS_c.txt -outb=$OUT_NAME-BARCODE_COUNTS_b.txt -outt=$OUT_NAME-TRANSCRIPT_COUNTS_t.txt

mv $HOME_DIR/$OUT_NAME-READ_COUNTS_c.txt $START_DIR/$OUT_NAME
mv $HOME_DIR/$OUT_NAME-BARCODE_COUNTS_b.txt $START_DIR/$OUT_NAME
mv $HOME_DIR/$OUT_NAME-TRANSCRIPT_COUNTS_t.txt $START_DIR/$OUT_NAME
mv $START_DIR/$OUT_NAME.gene.cout.csv $START_DIR/$OUT_NAME/$OUT_NAME.gene.cout.csv

/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME/$OUT_NAME.sam > $START_DIR/$OUT_NAME/$OUT_NAME.bam

# Remove unimportant files
rm $START_DIR/$OUT_NAME/$OUT_NAME.sam
rm $START_DIR/$LOG.sai
rm $START_DIR/$LOG2.sai
rm $START_DIR/$FASTQ_R1.all
rm $START_DIR/$FASTQ_R2.all
rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2


#copy all important files to a single location for that run
TRANSCRIPTS="_TranscriptFiles"
COUNTS="_CountFiles"
BARCODE="_BarcodeFiles"
SOUT="_soutFiles"

mkdir $PAST_DIR/$RUN_NAME$TRANSCRIPTS
cp $START_DIR/$OUT_NAME/$OUT_NAME-TRANSCRIPT_COUNTS_t.txt $PAST_DIR/$RUN_NAME$TRANSCRIPTS

mkdir $PAST_DIR/$RUN_NAME$COUNTS
cp $START_DIR/$OUT_NAME/$OUT_NAME-READ_COUNTS_c.txt $PAST_DIR/$RUN_NAME$COUNTS

mkdir $PAST_DIR/$RUN_NAME$BARCODE
cp $START_DIR/$OUT_NAME/$OUT_NAME-BARCODE_COUNTS_b.txt $PAST_DIR/$RUN_NAME$BARCODE

mkdir $PAST_DIR/$RUN_NAME$SOUT
cp $START_DIR/$OUT_NAME/$OUT_NAME.sout $PAST_DIR/$RUN_NAME$SOUT