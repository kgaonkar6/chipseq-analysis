#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


my($input, $sample, $filter, $help,$outvcf,$outrsID);
GetOptions("in|i=s"	=> \$input,
	"filter|f=s"	=> \$filter,
	"outvcf|o=s"	=> \$outvcf,
	"outrsID|r=s"	=> \$outrsID,
	"drugcatalog|d=s"=> \$drugcatalog,
	"help|h|?"	=> \&help);


open FH, "$input" or die " opening $input : $! \n";
print $outvcf."\t".$outrsID;
open OUT_VCF, ">","$outvcf" or die "opening file: $outvcf";
open OUT_rsID, ">","$outrsID" or die "opening file: $outrsID";

my @cols;
my @a;
my $rsID_col;

while (my $l = <FH>)	{
	chomp $l;
	my $filter_pass=0;
	if ($l =~ /#CHROM/) {
		print OUT_VCF $l."\n";
		@a = split(/\t/,$l);
		for (my $i =0; $i <scalar(@a);$i++){
		if ($a[$i] =~ /DP/){
		#print $a[$i];
		#print $i;
		push (@cols,$i);}
		if($a[$i] =~/dbSNP_142_GRCh38.ID/)
		{$rsID_col=$i;
		}
		}
		}
	else{ 
		@a = split(/\t/,$l);	
		for (my $i =0; $i <scalar(@cols);$i++){
		if(!($a[$cols[$i]]=~ /\./) && ($a[$cols[$i]] >= 10))
		{
		$filter_pass++;
		}
		}
		if ($filter_pass ==scalar(@cols))
		{
		print OUT_VCF $l."\n";
		if( !($a[$rsID_col] =~ /\./)){
		print OUT_rsID $a[$rsID_col]."\n";	}
		}
				
	}
}

close FH;		
close OUT_VCF;
close OUT_rsID;

`for ID in $(cat $outrsID ); do grep -w $ID variants.pharmgkb.tsv|cut -f4|cut -d "," -f1|cut -d ":" -f2; done >${outrsID}.filtered`
open FILT_rsID ,${outrsID}.filtered, or die "opening ${outrsID}.filtered: $!\n"
open FILT_rsIDJOB, ${outrsID}.filtered.job,or die "opening ${outrsID}.filtered.job: $!\n" 

while (my $ID =<FILT_rsID> ){
print FILT_rsIDJOB "echo". $ID ."|bior_lookup -d /data5/bsi/catalogs/user/v1/pharmacogenetics/pharmgkb/2013_07_03/variants.pharmgkb.tsv.bgz -p RSID |bior_drill -p DRUGS -p DISEASES";
}
close FILT_rsIDJOB;

`paste <(for id in $(source ${outrsID}.filtered.job |grep -v "#" |cut -f1);do awk -v ID=$id '{if( $21 == ID) print $0}' test.xls ;done ) <(source ${outrsID}.filtered.job |grep -v "#" |cut -f2-) > ${outrsID}.filtered.xls`




