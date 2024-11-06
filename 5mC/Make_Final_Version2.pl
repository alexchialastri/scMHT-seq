# Make a txt file that is compatible to earlier versions from AvO's *.faba file
# AvO's script has a bug in the exact nucleotide positions. This is corrected in this script

#***FOR MOUSE GENOME***
use Scalar::Util qw(looks_like_number);

$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCG      = substr($filename,0,-5)."_Full5mC_CG_Final_Rmdup.txt";
open($outputCG, '>', $outputfilenameCG);
$outputfilenameCG_sort = substr($filename,0,-5)."_Full5mC_CG_Final_Rmdup_sort.txt";

$outputfilenameCA      = substr($filename,0,-5)."_Full5mC_CA_Final_Rmdup.txt";
open($outputCA, '>', $outputfilenameCA);
$outputfilenameCA_sort = substr($filename,0,-5)."_Full5mC_CA_Final_Rmdup_sort.txt";

$outputfilenameCC      = substr($filename,0,-5)."_Full5mC_CC_Final_Rmdup.txt";
open($outputCC, '>', $outputfilenameCC);
$outputfilenameCC_sort = substr($filename,0,-5)."_Full5mC_CC_Final_Rmdup_sort.txt";

$outputfilenameCT      = substr($filename,0,-5)."_Full5mC_CT_Final_Rmdup.txt";
open($outputCT, '>', $outputfilenameCT);
$outputfilenameCT_sort = substr($filename,0,-5)."_Full5mC_CT_Final_Rmdup_sort.txt";

#open($file, "/hpc/hub_oudenaarden/s.dey/5mC_RNA_Emb_Hyb_CAST_B6_sc_20161116/HybEmb5mCRNA-Lib1_5mCReads/AlleleSpecificData/FinalManuscript/HybEmb5mCRNA-Lib1_5mCReads_B6_se.allele.faba");
#open($output, ">/hpc/hub_oudenaarden/s.dey/5mC_RNA_Emb_Hyb_CAST_B6_sc_20161116/HybEmb5mCRNA-Lib1_5mCReads/AlleleSpecificData/FinalManuscript/HybEmb5mCRNA-Lib1_5mCReads_B6_se_Full5mC_CG_Final_Rmdup.txt");

while ($line = <$file>)
{
    @x = split(" ",$line);

if (looks_like_number($x[1])) {
    $x[0] = $x[0] + 1;
    $x[1] = $x[1] + 1;
    if ($x[3] == 1)
    {
	$x[2] = $x[2] + 1;
    }
    if ($x[3] == -1)
    {
        $x[2] =$x[2] +3;
    }
    
    $y = substr($x[5],0,2);
    if ($y eq 'CG')
    {
	#print $outputCG ("$x[0] $x[1] $x[2] $x[4] $x[3] $x[5]\n");
	print $outputCG ("$x[0] $x[1] $x[2] $x[4] $x[3]\n");
    }
    
    if ($y eq 'CA')
    {
        #print $outputCA ("$x[0] $x[1] $x[2] $x[4] $x[3] $x[5]\n");
	print $outputCA ("$x[0] $x[1] $x[2] $x[4] $x[3]\n");
    }
    
    if ($y eq 'CC')
    {
        #print $outputCC ("$x[0] $x[1] $x[2] $x[4] $x[3] $x[5]\n");
	print $outputCC ("$x[0] $x[1] $x[2] $x[4] $x[3]\n");
    }
    
    if ($y eq 'CT')
    {
        #print $outputCT ("$x[0] $x[1] $x[2] $x[4] $x[3] $x[5]\n");
	print $outputCT ("$x[0] $x[1] $x[2] $x[4] $x[3]\n");
    }
}
}

system("sort -k1,1n -k2,2n -k3,3n $outputfilenameCG > $outputfilenameCG_sort");
 
system("sort -k1,1n -k2,2n -k3,3n $outputfilenameCA > $outputfilenameCA_sort");

system("sort -k1,1n -k2,2n -k3,3n $outputfilenameCC > $outputfilenameCC_sort");

system("sort -k1,1n -k2,2n -k3,3n $outputfilenameCT > $outputfilenameCT_sort");
