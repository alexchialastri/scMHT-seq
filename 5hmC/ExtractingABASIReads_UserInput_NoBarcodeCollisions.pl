# Extract 5hmC reads based on 96 barcodes from fastq files
                                                                                                                
use warnings;
use Getopt::Long;
 
GetOptions (   
	    "FASTQ_R1=s" => \$FASTQ_R1,
        "CELSEQ_BC=s" => \$BARCODES,
		"scMHseq=i" => \$scMHseq,
		"MSPJI_BC=s" => \$MSPJIBCpath,
		"scTHseq=i" => \$scTHseq,
		"ABASI_BC=s" => \$ABASIBCpath
    )   
 or die("Error in command line arguments\n usage: All Arguments in this order are needed: --FASTQ_R1 is the input name of your input R1 fastq file && --CELSEQ_BC is a .csv file of Cel-seq cell specific barcodes RNA && --scMHseq is a true (1) flase (0) interger to search for barcode collisions between ABASI and MSPJI sequencing && MSPJI_BC is the file path to the MSPJI barcodes && --scTHseq is a true (1) flase (0) interger to search for barcode collisions between RNA and ABASI sequencing && ABASI_BC is the file path to the ABASI barcodes )\n");
 
 
open($fastafileR1, "$FASTQ_R1");

my $outputR1In = substr($FASTQ_R1,0,-6)."-ABA.fastq";

open($outputR1, '>', "$outputR1In");


# READ in MSPJI Barcodes
if ($scMHseq == 1) {
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



# Read in CEL SEQ Barcodes if is part of joint sequencing
if ($scTHseq == 1) {
	open($BC, "$BARCODES"); #opens celseq barcode file
	$i = 0;
		while ($BCline = <$BC>)
		{
			chomp $BCline;   
			@BarC = split("\t",$BCline);
			$BarCode[$i] = $BarC[1];
			$i++;     
		}
}
    

# Read in ABASI Barcodes
open($ABASIBCfile, "$ABASIBCpath"); #opens ABASI barcode file
$i = 0;
while ($BCline3 = <$ABASIBCfile>)
 {
     chomp $BCline3;   
     @ABASIBarC = split(";",$BCline3);
     $ABASIBarCode[$i] = $ABASIBarC[0];
     $i++; 
 }
    


while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;


    $RNABCSection = substr($sequence,0,8); #Get CELSEQ Barcode Position sequence
    $MSPJIBCSection = substr($sequence,3,8); #Get MSPJI Barcode Position sequence
    $ABASIBCSection = substr($sequence,0,6); #Get ABASI Barcode Position sequence
	
    
	
	
	$RNACount = 0;
	$MSPJICount = 0;
	$ABASICount = 0;

        #Check if ABASI Barcode is a match
	foreach (@ABASIBarCode){
		if ($_ eq $ABASIBCSection){$ABASICount++;}
       	}	
	
		
	
	
	if ($scTHseq == 1) {
	$PosPolyT = substr($sequence,12,5); #Tie Breaker to check if there is a 5 bp polyT sequence, which should be true for CELSEQ
		foreach (@BarCode){ #Check if CELSEQ Barcode is a match
			if ($_ eq $RNABCSection){$RNACount++;}
		}
	}
	else {$PosPolyT = "GGGGG";} #unmatching PosPolyT seq if not combined with RNA }
	

	
	if ($scMHseq == 1) { 
	    #Check if MSPJI Barcode is a match
	    foreach (@MSPJIBarCode){
		if ($_ eq $MSPJIBCSection){$MSPJICount++;}
	    }
	    
	}
	
	$NonABASISEQ = $RNACount + $MSPJICount; #Match count of all non AbaSI Barcodes
	
   # print "$MSPJICount $NonABASISEQ $ABASICount\n";
	#Print out New Fastq File if barcode match only AbaSI
	if ( $ABASICount > 0){
		if ( $NonABASISEQ == 0 || ($MSPJICount == 0 && $PosPolyT ne "TTTTT")){
			print $outputR1 ("$header$sequence$plus$qual");
		#	print "output printed\n";
		}
	}
	

}
