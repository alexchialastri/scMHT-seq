#!/usr/bin/perl -w -s
use lib '/home/achialastri/perlscripts/mRNA_Mapping';
use tools;

if (!($r && $f1 && $out && $t)){
    die "usage: -r=REFERENCE -f1=READ1 -f2=READ2(optional) -out=OUTPUT_PREFIX -outdir=OUTPUT_DIRECTORY (optional) -t=THREADS -ind=is or bwtsw (default: is (<2GB reference) -q=PHRED_QUALITY_FOR_TRIMMING (0..62) -aln_n=edit distance -aln_k=edit distance in seed -l=SEED_LENGTH -BL=barcode length left read (bwa aln option -B, default is 0) -BR=barcode length right read (bwa aln option -B, default is 0) -i= 1 or 0 (1 if indexing is required) -gff=DATA.gff (optional, produces wig files) -s_flag=1 (optional: strand specific mapping) -u=1 (optional: only reads mapping to one strand) -uniq=1 (optional, keep only uniquely mapped reads) -test=1 -npr= 0 or 1 (1: do not process_sam.pl 2: run only do not process_sam.pl ) -nsam= 0 or 1 (1: do not produce new sam file) -cel=0 or 1 (1: process CEL-seq *.sam output) -bar=cel-seq_barcodes.csv -fstr= 1 or 0 ( if 1 only mappings to the sense strand are allowed ) -anno= anno.csv ( optional (CEL-seq), when mapping to the genome, run get_anno.pl on sam file first) -rb= 0 or 1 (optional (CEL-seq), use random barcodes) -rb_len= length of random barcode (default = 4) -dprm= 0 or 1 (for CEL-seq: 1: reomve pcr duplicates)\n";
}
# hard coded:
$n = 100; # maximal number of unpaired reads to output in XA tag
$N = 100; # maximal number of paired reads to output in XA tag
$q = 0 if !$q; # base quality cutoff for trimming
$l = 200 if !$l; # seed length
$i = 0 if !$i; # build reference database
$BL = 0 if !$BL;
$BR = 0 if !$BR;
$npr  = 0 if !$npr;
$nsam = 0 if !$nsam;
$cel = 0 if !$cel;
$fstr = 0 if !$fstr;
$ind = "is" if !$ind;
$anno = 0 if !$anno;
$rb = 0 if !$rb;
$rb_len = 4 if !$rb_len;
$dprm = 0 if !$dprm;

$aln_n = 0.04 if !$aln_n; # edit distance
$aln_k = 2 if !$aln_k; # edit distance in seed

$test = 0 if !$test;
if ($outdir){
    makedir($outdir);
    $out = $outdir."/".$out;
   
}
 print $out;
$pflag = 0;
$pflag = 1 if ($f1 && $f2);


if ($i){
    $str = "/home/cwangsanuwat/bwa/bwa-0.7.15/bwa index -a ".$ind." ".$r;
    print $str."\n";
    (system($str) == 0 or die "Could not execute ".$str."\n") if ($test == 0);
}

@F = ($f1);
@F = (@F, $f2) if $pflag;

@G = @F;
@H = @F;

for $i (0..$#F){
    $G[$i] =~ s/(\.)\w+$/\.fastq/;
    $H[$i] =~ s/(\.)\w+$/\.sai/;
}

if ( $npr != 2 ){
    for $i (0..$#G){
	if ($F[$i] =~ /txt/){
	    $str = "qseq2fastq.pl -clean=1 -in=".$F[$i]." > ".$G[$i];
	    print $str."\n";
	    (system($str) == 0 or die "Could not execute ".$str."\n") if ($test == 0);
	}
	if ( $i == 0 ) { $B = $BL; } else { $B = $BR; }
	$str = "/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -B ".$B." -q ".$q." -n ".$aln_n." -k ".$aln_k." -l ".$l." -t ".$t." ".$r." ".$G[$i]." > ".$H[$i];
	print $str."\n";
	(system($str) == 0 or die "Could not execute ".$str."\n") if ($test == 0);
    }
    
    if ( $nsam == 0 ){
	if ($pflag){
	    $str = "/home/cwangsanuwat/bwa/bwa-0.7.15/bwa sampe -n ".$n." -N ".$N." ".$r." ".join(" ",(@H, @G))." > ".$out.".sam";
	}else{
	    $str = "/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n ".$n." ".$r." ".join(" ",(@H, @G))." > ".$out.".sam";
	}
    }
    print $str."\n";
    (system($str) == 0 or die "Could not execute ".$str."\n") if ($test == 0);
    
}
if ( $npr == 0 || $npr == 2){
    $s = 1;
    $s = 0 if $pflag;
    $str = "/home/achialastri/perlscripts/mRNA_Mapping/process_sam_strand_Deylab.pl -in=".$out.".sam -s=".$s;
    if ( $cel ){
      $str = "/home/achialastri/perlscripts/mRNA_Mapping/process_sam_cel_seq_Deylab.pl -in=".$out.".sam -bc=".$bar." -anno=".$anno." -rb=".$rb." -rb_len=".$rb_len;
      if ( $fstr ){
	$str = $str." -fstr=".$fstr;
      }
      if ( $dprm ){
	$str = $str." -dprm=".$dprm;
      }
    }
    if ($gff){
	$str = $str." -gff=".$gff;
    }
    if ($s_flag){
	$str = $str." -s_flag=".$s_flag;
    }
    if ($u){
	$str = $str." -u=".$u;
    }
    if ($uniq){
	$str = $str." -uniq=".$uniq;
    }
    print $str."\n";
    (system($str) == 0 or die "Could not execute ".$str."\n") if ($test == 0);
}
