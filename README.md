# scMHT-seq
Analysis Scripts for scMHT-seq

scMHT-seq contains two sequencing libraries, one library that contains all molecule types which is dominated by gDNA derived molecules. This fastq file is used for 5mC and 5hmC analysis. Prior to running the shell script for 5mC or 5hmC (found in their own folders), read 1 of this fastq file is trimmed to 76 basepairs using the MakeFastqShorter.pl script. The other sequencing library is derived from an RNA enriched library. Prior to running the shell script for RNA (found in the RNA folder), read 1 of this fastq file is trimmed to 25 basepairs using the MakeFastqShorter.pl script and read 2 was trimmed to 50 basepairs using the MakeFastqShorter.pl script. These trimmed fastq files are the input into the pipelines. The shell script pipelines need to be updated to contain the directories containing the perl scripts and other tools used in the pipeline (clumpify from bbmap, bwa, samtools..etc.).



