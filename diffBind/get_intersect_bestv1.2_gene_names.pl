#!/usr/bin/perl
use strict;
use warnings;
## Perl script to get uniq geneID@gene list from per transcript analysis
## Krutika Gaonkar 1/27/2016
use Data::Dumper;

if (scalar(@ARGV) != 4)	
{
	die ( "USAGE: get_intersect_best.pl [transcript] [DE.list] [state] [output dir]\n" );
}


my $transcript = $ARGV[0];
my $DE_list= $ARGV[1];
my $state=$ARGV[2];
my $dir=$ARGV[3];


my $sample=`basename $transcript|cut -d "." -f1`;
chomp ($sample);

open BED, "$dir/$transcript" or die;
open BEST, ">" ,"$dir/final_max$sample" or die;

my %cons_DE_pertrans_info = ();
our %cons_DE_info=();
### chr_merged.start_merged.end -> [best FDR/qvalue]


open DE ,"$DE_list" or die;

while (<DE>){
        my $row = $_;
        chomp $row;
        our @line = split("\t",$row);
	#my @gene_name=split("@",$line[1]);
		#print $line[4];
 	       if ( !$cons_DE_info{"$line[3]"}) {$cons_DE_info{"$line[3]"}="$line[4]_$line[5]";}
	#	print $gene_name[0];
}

while (<BED>){
 	       	my $row = $_;
       		chomp $row;
	        my @line2 = split("\t",$row);

#		$cons_DE_pertrans_info{"$line2[3]"}="$line2[6]_$line2[7]_$line2[8]_$line2[9]_$line2[10]_$line2[11]";
		$cons_DE_pertrans_info{"$line2[0]"}="$line2[2]_$line2[3]_$line2[4]_$line2[5]_$line2[6]_$line2[7]";
		#print $line2[0];
}



foreach my $geneDE ( sort keys %cons_DE_info)
{
our $value="NA_NA_NA_NA_NA_NA";
our $value_2=$cons_DE_info{$geneDE}; 
$value_2=~s/_/\t/g;
foreach my $geneID (sort keys %cons_DE_pertrans_info)
{
my @gene=split("@",$geneDE);
if ($geneID=~/^$gene[0]/)
{
$value=$cons_DE_pertrans_info{$geneID};
$value =~s/_/\t/g;

print BEST "$geneID\t$value\t$value_2\n";
}
}
}
###`rm $dir/intersect_$name`;

