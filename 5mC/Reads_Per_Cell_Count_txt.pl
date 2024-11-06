# Reads per cell from a txt file that was derived from FABA file (aka CG, CA, CC, CT files)



$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCG      = substr($filename,0,-4)."_CountsPerCell.txt";
open($outputCG, '>', $outputfilenameCG);


for ($j=0;$j<=95;$j++)
    {
	$ReadsPerCell[$j] = 0;

    }
	
while ($line = <$file>)
{
    @x = split("\t",$line);
    $ReadsPerCell[$x[0]-1]++;
	#$ReadsPerCell[$x[0]]++;
	
}

for ($j=0;$j<=95;$j++)
    {
print $outputCG ("$ReadsPerCell[$j]\n");
}
