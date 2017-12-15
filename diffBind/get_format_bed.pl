#!/usr/bin/perl
use strict;
use warnings;
## Perl script to generate the files required for diffBind
## Krutika Gaonkar 1/27/2016
use Data::Dumper;

if (scalar(@ARGV) != 3)	
{
	die ( "USAGE: get_intersect_best.pl [bed for mark] [caller] [output dir]\n" );
}


my $peak_bed = $ARGV[0];
my $caller= $ARGV[1];
my $dir=$ARGV[2];


my $sample=`basename $peak_bed|cut -d "." -f1`;
chomp ($sample);
my $name=$sample."_formated.bed";
open BED, "$peak_bed" or die;
open FORMAT_BED, ">" , "$dir/$name" or die;

my %cons_peak_info = ();
### chr_merged.start_merged.end -> [-log10(FDR)/-log10(qvalue)]

my $FDR=0;


if ( $caller eq "sicer" )
{
while (<BED>){
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
	if ($line[7] == 0){$FDR=0;}
	else{
	$FDR=-log10($line[7]);
	}
	#print $raw_FDR."\t".$FDR."\t";	
        $cons_peak_info{"$line[0]_$line[1]_$line[2]"} = "$FDR";
}
}

if ( $caller eq "macs" )
{
while (<BED>){
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
	my $q_val=$line[8];
        $cons_peak_info{"$line[0]_$line[1]_$line[2]"} = "$q_val";
}
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}



foreach my $region ( sort keys %cons_peak_info)
{
my $value=$cons_peak_info{$region};
$region =~s/_/\t/g;
print FORMAT_BED "$region\t$value\n";
}

###`rm $dir/intersect_$name`;

