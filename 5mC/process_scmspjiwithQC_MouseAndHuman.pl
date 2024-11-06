#!/usr/bin/perl

use strict "subs";
use strict "refs";
#use strict "vars";

print STDERR "
  processes sam file of single-cell MspJI-seq data
 
  Usage
  process_scaba.pl genome_file.fa sam_file.sam msjp1_barcodes.csv
 
";

# read in genome file
my $genomefile = $ARGV[0];
open(MULTIFASTA, "$genomefile") || die " FATAL ERROR:\n Unable to load '$genomefile'.\n";
while (<MULTIFASTA>) {
	chomp($_);
	if ($_=~/^>(.*)/) {
		if ($seq) {
			$sequence{$header}=$seq;
		}
		$header = $1;
		$counter++;
		$seq    = '';
	} else {
		$seq.=$_;
	}
 
}
$sequence{$header}=$seq;

# read in cell-specific barcode file
my $file2 = $ARGV[2] or die "Need to get cell specific barcode file on the commend line\n";
my @csbc;
open(my $cscodes, '<', $file2) or die "Could not open '$file2' $!\n";
while (my $line = <$cscodes>) {
  chomp $line;
    my @fields = split "\t" , $line;
    push @csbc,@fields[0];
}

# read in sam file, filter reads, and count CG dinucleotides
my $file1 = $ARGV[1] or die "Need to get SAM file on the command line\n";
open(my $data, '<', $file1) or die "Could not open '$file1' $!\n";

# for QC file
for ($k=0;$k<96;$k++)
 {
   $allreads[$k]=0;
   $rawreads[$k] = 0;
   $mappedreads[$k] = 0;
   $cleanreads[$k] = 0;
 }


my $rABAoutputfile = substr($file1,0,-3)."raba";
my $fABAoutputfile = substr($file1,0,-3)."faba";
my @chrvec=qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM ZymoDNAStandardSet PhageLambda pUC19 MM10_chr1 MM10_chr2 MM10_chr3 MM10_chr4 MM10_chr5 MM10_chr6 MM10_chr7 MM10_chr8 MM10_chr9 MM10_chr10 MM10_chr11 MM10_chr12 MM10_chr13 MM10_chr14 MM10_chr15 MM10_chr16 MM10_chr17 MM10_chr18 MM10_chr19 MM10_chrX MM10_chrY MM10_chrM);

open(my $abaout, '>', $rABAoutputfile) or die "Can't open file for writing: $!\n";

$N_plusmatch=0;
$N_minmatch=0;

while (my $line = <$data>) {
    chomp $line;
    my @fields = split "\t" , $line;
    my $seqflag = $fields[1];
    my $chr = $fields[2];
    my $cutcoord = $fields[3];
    my $cigar = $fields[5];
    my $cellbarcode = substr($fields[11],8,8);
    my $umi = substr($fields[11],5,3);
    my $uniqueflag = $fields[12];
    my $strand = 0;
    my $offset = 0;
    my $CGcount = 0;
 
    if ($cellbarcode ~~ @csbc ){

	$allreads[$csbcindex]++;}

    if ($cellbarcode ~~ @csbc & $chr ~~ @chrvec){
	my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
	$rawreads[$csbcindex]++;}
    if (($seqflag == 0 | $seqflag == 16) & ($cigar =~ /(^\d\d)M/) & (length($cigar)== 3) & $uniqueflag =~ "XT:A:U" & $cellbarcode ~~ @csbc & $chr ~~ @chrvec){
       my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
       my ($chrindex) = grep {$chrvec[$_] ~~ $chr} 0 .. 49;       
       $mappedreads[$csbcindex]++;
       my $cutsite = uc(substr($sequence{$chr},$cutcoord-20,125));
       my $Nhits=0;

       my $MatchLength = substr($cigar,0,2);
       my $AdjustedLength = 65 - int($MatchLength);
      
		#Indirect, 5mC originally on + Strand 
         if ($seqflag == 0 & substr($cutsite,6,1) eq 'C'){
		 $CGpos=$cutcoord-14; $CGcand = uc(substr($sequence{$chr},$CGpos,4));$strand=1; $Nhits++;
		 }
		 
		 #Direct, 5mC originally On + strand
         if ($seqflag == 16 & substr($cutsite,67-$AdjustedLength,1) eq 'C'){
		 $CGpos=$cutcoord+47-$AdjustedLength; $CGcand = uc(substr($sequence{$chr},$CGpos,4));$strand=1; $Nhits++;
		 }
         
		 #Direct, 5mC originally on - strand
		 if ($seqflag == 0 & substr($cutsite,35,1) eq 'G'){
		 $CGpos=$cutcoord+12; $CGcand = uc(substr($sequence{$chr},$CGpos,4)); $seqmatch = reverse $CGcand; $seqmatch =~ tr/ACGTacgt/TGCAtgca/; $CGcand = $seqmatch; $strand=-1; $Nhits++;
		 }
		
		 #Indirect, 5mC orignally on - strand 
		 if ($seqflag == 16 & substr($cutsite,96-$AdjustedLength,1) eq 'G'){
		 $CGpos=$cutcoord+73-$AdjustedLength; $CGcand = uc(substr($sequence{$chr},$CGpos,4));$seqmatch = reverse $CGcand; $seqmatch =~ tr/ACGTacgt/TGCAtgca/; $CGcand = $seqmatch; $strand=-1; $Nhits++;
		 }
#note for last two 'G' CGpos is shifted by 3bp because of reverse complement !!    




       if ($Nhits == 1){
	   if ($chrindex == 25) {$chrindex="Zymo";}
	   if ($chrindex == 26) {$chrindex="PhageLambda";}
	   if ($chrindex == 27) {$chrindex="pUC19";}
	   if ($chrindex > 27) {$tempChrindex = $chrindex-27; $chrindex="MM10_".$tempChrindex;}
	   print $abaout "$csbcindex\t$chrindex\t$CGpos\t$strand\t$umi\t$CGcand\t$Nhits\t$seqflag\n";
$cleanreads[$csbcindex]++;};
}} 

# write QC count to file
my $QCoutputfile = substr($file1,0,-3)."QC";
open(my $qcout, '>', $QCoutputfile) or die "Can't open file for writing: $!\n";
my $QCoutputfile2 = substr($file1,0,-3)."QC2";
open(my $qcout2, '>', $QCoutputfile2) or die "Can't open file for writing: $!\n";
for ($k=0;$k<96;$k++)
 {
   print $qcout "$k \t $rawreads[$k] \t $mappedreads[$k] \t $cleanreads[$k] \n";
   print $qcout2 "$allreads[$k]\n";
 }


# keep only unique cuts and sort
system("sort -u $rABAoutputfile > $fABAoutputfile");

