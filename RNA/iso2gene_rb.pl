#!/usr/bin/perl -w -s

#use lib "/Users/d.grun/data/bin";
use lib "/home/achialastri/perlscripts/mRNA_Mapping";
use tools;

if (scalar @ARGV == 1){
    die "usage: -m=iso2gene.csv -in=iso_counts.csv (comma separated list) -out=aggr_counts.csv\n" if $ARGV[0] eq "help";
}
open(IN,"<",$m);
while(<IN>){
  chomp;
  ($g,$i) = split(/\t/);
  push(@{$h{$g}}, split(/\|/,$i) );
}
close(IN);

@ina = split(/\,/,$in);

foreach $k (@ina){
  open(IN,"<",$k);
  while(<IN>){
    if ( $_ =~ /GENEID/ ){
      $title = $_;
      next;
    }
    chomp;
    @F = split(/\t/);
    if ( !exists($r{$F[0]}{$F[1]}) ){
      push(@{$r{$F[0]}{$F[1]}}, @F[2..$#F] );
    }else{
      for $j (2..$#F){
	${$r{$F[0]}{$F[1]}}[$j - 2] += $F[$j];
      }
    }
  }
  close(IN);
}



open(OUT,">",$out);
print OUT $title;
foreach (sort keys %h){
  %G = ();
  foreach $i (@{$h{$_}}){
    next if !exists($r{$i});
    foreach $b ( sort keys %{$r{$i}} ){
      $seen{$i} = 1;
      if ( !exists($G{$b}) ){
	@{$G{$b}} = @{$r{$i}{$b}}
      }else{
	for $k ( 0..$#{$G{$b}}){
	  ${$G{$b}}[$k] = ${$G{$b}}[$k] + ${$r{$i}{$b}}[$k];
	}
      }
    }
  }
  foreach $b (sort keys %G){
    print OUT join("\t",($_,$b,@{$G{$b}}))."\n";
  }
}
foreach (sort keys %r){
  foreach $b ( sort keys %{$r{$_}} ){
    print OUT join("\t",($_,$b,@{$r{$_}{$b}}))."\n" if !exists($seen{$_});
  }
}
close(OUT);
