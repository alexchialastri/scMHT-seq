# Shorten a fastq file i.e. 150 read length to our normal 76
                                                                                                                
use warnings;

$Input = $ARGV[0];
$OutputName = $ARGV[1];
$LengthMake = $ARGV[2]; 

if (!($ARGV[0] && $ARGV[1] && $ARGV[2])){
    die "usage: All Arguments in this order are needed: Argument 0 is the input fastq file && Argument 1 is the desired trimmed output fastq file name && Argument 2 is the number of base pairs to trim to\n";
}

open($fastafileR1, $Input);
open($outputR1,'>',$OutputName);


while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;

    $sequence2 = substr($sequence,0,$LengthMake);
    $qual2 = substr($qual,0,$LengthMake);
    
    print $outputR1 ("$header$sequence2\n$plus$qual2\n");
	    
	
    
}
