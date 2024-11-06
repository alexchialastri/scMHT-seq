# Extract CelSeq reads based on 96 barcodes from fastq files (containing CEL-Seq reads)                                                                                                 
use warnings;
use Getopt::Long;
 
GetOptions (   
              "StartingDir=s"   => \$START_DIR,      # string
	    "OUTNAME_R1=s" => \$OUT_NAME_R1,
	     "OUTNAME_R2=s" => \$OUT_NAME_R2,
	    "FASTQ_R1=s" => \$FASTQ_R1,
	    "FASTQ_R2=s" => \$FASTQ_R2,
            "CELSEQ_BC=s" => \$BARCODES,
    "scMATseq=i" => \$scMATseq,
    "MSPJI_BC=s" => \$MSPJIBCpath,
    "scTHseq=i" => \$scTHseq,
    "ABASI_BC=s" => \$ABASIBCpath
    )   
 or die("Error in command line arguments\n usage: All Arguments in this order are needed: --StartingDir is the starting directory && --OUTNAME_R1 is the output name of your output R1 fastq file && --OUTNAME_R2 is the output name of your output R2 fastq file && --FASTQ_R1 is the input name of your input R1 fastq file && --FASTQ_R2 is the input name of your input R2 fastq file && --CELSEQ_BC is a .csv file of Cel-seq cell specific barcodes RNA && --scMATseq is a true (1) flase (0) interger to search for barcode collisions between CELSEQ and MSPJI sequencing && MSPJI_BC is the file path to the MSPJI barcodes && --scTHseq is a true (1) flase (0) interger to search for barcode collisions between CELSEQ and ABASI sequencing && ABASI_BC is the file path to the ABASI barcodes )\n");



#Get full file names and make new ones for output
$FileName1=join "",$START_DIR,"/",$FASTQ_R1;
$FileName2=join "",$START_DIR,"/",$FASTQ_R2;
$FileName3=join "",$START_DIR,"/",$OUT_NAME_R1,".fastq";
$FileName4=join "",$START_DIR,"/",$OUT_NAME_R2,".fastq";



open($fastafileR1, "$FileName1"); #Opens Fastq file
open($fastafileR2, "$FileName2");

open($outputR1, ">$FileName3"); #opens to write into output file
open($outputR2, ">$FileName4");

open($BC, "$BARCODES"); #opens celseq barcode file




# Read in CEL SEQ Barcodes
$i = 0;
while ($BCline = <$BC>)
{
    chomp $BCline;   
    @BarC = split("\t",$BCline);
    $BarCode[$i] = $BarC[1];
    $i++;     
}

# READ in MSPJI Barcodes if is part of joint sequencing
if ($scMATseq == 1) {
    open($MSPJIBCfile, "$MSPJIBCpath"); #opens MSPJI barcode file
    $i = 0;
    while ($BCline2 = <$MSPJIBCfile>)
    {
    chomp $BCline2;   
    @MSPJIBarC = split("\t",$BCline2);
    $MSPJIBarCode[$i] = $MSPJIBarC[0];
    $i++;
    }
 }       

# Read in ABASI Barcodes if is part of Joint sequencing
if ($scTHseq == 1) {
    open($ABASIBCfile, "$ABASIBCpath"); #opens ABASI barcode file
    $i = 0;
    while ($BCline3 = <$ABASIBCfile>)
    {
    chomp $BCline3;   
    @ABASIBarC = split(";",$BCline3);
    $ABASIBarCode[$i] = $ABASIBarC[0];
    $i++; 
    }
 }                        

while ($header = <$fastafileR1>)
{
#Read in all parts of the 4 line FASTQ format from both R1 and R2
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;

    $headerR2   = <$fastafileR2>;
    $sequenceR2 = <$fastafileR2>;
    $plusR2     = <$fastafileR2>;
    $qualR2     = <$fastafileR2>;

    $RNABCSection = substr($sequence,0,8); #Get CELSEQ Barcode Position sequence
    $MSPJIBCSection = substr($sequence,3,8); #Get MSPJI Barcode Position sequence
    $ABASIBCSection = substr($sequence,0,6); #Get ABASI Barcode Position sequence
    $PosPolyT = substr($sequence,14,5); #Tie Breaker to check if there is a 5 bp polyT sequence, which should be true for CELSEQ. Made 14 to account for 6 UMI reads we have now 2/24/21
	
	
	$RNACount = 0;
	$MSPJICount = 0;
	$ABASICount = 0;
	
	foreach (@BarCode){ #Check if CELSEQ Barcode is a match
		if ($_ eq $RNABCSection){$RNACount++;}
		}
	
	
	if ($scMATseq == 1) { #Check if MSPJI Barcode is a match
	foreach (@MSPJIBarCode){
		if ($_ eq $MSPJIBCSection){$MSPJICount++;}
		}
	}
	
	if ($scTHseq == 1) { #Check if ABASI Barcode is a match
	foreach (@ABASIBarCode){
		if ($_ eq $ABASIBCSection){$ABASICount++;}
		}
	}
	
	$NonCELSEQ = $MSPJICount + $ABASICount;
	
  
	if ( $RNACount > 0){
		#print "$PosPolyT\t$NonCELSEQ\n";
		if ( $NonCELSEQ == 0 || $PosPolyT eq "TTTTT"){
			print $outputR1 ("$header$sequence$plus$qual");
			print $outputR2 ("$headerR2$sequenceR2$plusR2$qualR2");
		}
	}
   
}
