#!/usr/bin/perl


my %snp;
my $site;

open FH1, "$ARGV[0]" or die "";
open FH2, "$ARGV[1]" or die "";
open (my $FH3,'>>', $ARGV[2]) or die"";

$attribute= $ARGV[3];
print $attribute;

while ( $line = <FH1> )
{
	chomp($line);
	@info1=split(/\t/,$line);
	$key1=$info1[0].":".$info1[2];
	shift @info1;
	$snp{$key1}= "0";
}

while ( $line = <FH2> )
{
        chomp($line);
        @info2=split(/\t/,$line);
	if ( $info2[7] eq $attribute )
	{ 
	$site=$info2[0].":".$info2[2];	
	$snp{$site}="1";
	}
}

	
foreach my $key1 ( sort keys %snp)
{
print $FH3 "$key1\t$snp{$key1}\n";
}



