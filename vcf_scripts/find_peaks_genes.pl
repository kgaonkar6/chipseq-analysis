#!/usr/bin/perl
use strict;
use warnings;
## Perl script takes in peak positions and finds genes with nearby TSS or TTS and reports their distance
## Jared Evans 12/19/11
## modified Krutika Gaonkar 12/18/15
use Data::Dumper;

if (scalar(@ARGV) != 3)	
{
	die ( "USAGE: find_nearby_genes.pl [peaks BED file] [reflat file] [output file] [neighborhood distance] [min dist] \n" );
}

open VARIANTS, "$ARGV[0]" or die "opening file: $ARGV[0]";
open INFO, "$ARGV[1]" or die "opening file: $ARGV[1]";
open OUT, "+>","$ARGV[2]" or die "opening file: $ARGV[2]";





my %var = ();
## chr_start -> [peakcenter, peakrank, chr, start, stop, peaklength]
my %info = ();
## chr__start -> [INFO]

open OUT, "+>","$ARGV[0]" or die "opening file: $ARGV[0]";
# load peaks into hash
while (<VARIANTS>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $variant = $line[0]_$line[1];
	$var{$varinat} = ".\t.\t.\t.\t."
	#load all peaks once in output hash
}

# load transcripts into hash of hashes
while(<INFO>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $variants=$line[0]_$line[1];	
	$info{$line[0]_$line[1]} = [
}

print OUT "Peak_Center\tTranscript\tChr\tGene\tStart\tStop\tTSS\tTTS\tStrand\tINFO\n";

for (my $a=0; $a<scalar(@max_dist); $a++)
{
	
foreach my $peak_line (keys %peaks){
	#my $peak_line = "chr1_5011760_5014166"; #upstream -
	#my $peak_line = "chr1_15786175_15786773"; #upstream +
	#my $peak_line = "chr1_13355099_13355685"; #downstream -
	#my $peak_line = "chr1_34076808_34079752"; #downstream +
	my @peak_info = @{$peaks{$peak_line}};
	foreach my $transcript (keys %{$transcripts{$peak_info[2]}}){
		my @transcript_info = @{$transcripts{$peak_info[2]}{$transcript}};
		## TSS:
		if(($peak_info[0]-$transcript_info[4]) <= $max_dist[$a] && ($peak_info[0]-$transcript_info[4]) >= $min_dist[$a] && $transcript_info[4]>$transcript_info[5]){
				print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tUpstream:$max_dist[$a]"."\n";
				$type{"Peak_".$peak_info[0]."_Upstream:".$max_dist[$a]}=1;
		}elsif(($transcript_info[4]-$peak_info[0]) <= $max_dist[$a] && ($transcript_info[4]-$peak_info[0]) >= $min_dist[$a] && ($transcript_info[5]>$transcript_info[4])){
				print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tUpstream:$max_dist[$a]"."\n";
				$type{"Peak_".$peak_info[0]."_Upstream:".$max_dist[$a]}=1;
		}
		## TTS:
		if(($peak_info[0]-$transcript_info[5]) <= $max_dist[$a] && ($peak_info[0]-$transcript_info[5]) >= $min_dist[$a] && $transcript_info[5]>$transcript_info[4] ){
				print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tDownstream:$max_dist[$a]"."\n";
				$type{"Peak_".$peak_info[0]."_Downstream:".$max_dist[$a]}=1;
		}elsif(($transcript_info[5]-$peak_info[0]) <= $max_dist[$a] && ($transcript_info[5]-$peak_info[0]) >= $min_dist[$a] && $transcript_info[4]>$transcript_info[5]){
				print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tDownstream:$max_dist[$a]"."\n";
				$type{"Peak_".$peak_info[0]."_Downstream:".$max_dist[$a]}=1;
		}
		#UTR,exon and intron regions	
		if( $max_dist[$a] == 0  && $transcript_info[3] >= $peak_info[0] && $peak_info[0]>= $transcript_info[2]){	
		#5'UTR and 3'UTR
		if ( $transcript_info[6] eq "+" && $peak_info[0] >= $transcript_info[4] && $peak_info[0]<= $transcript_info[10] )
		{
		                print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\t5'UTR"."\n";
				$type{"Peak_".$peak_info[0]."_5'UTR"}=1;
		}elsif ( $transcript_info[6] eq "-" && $peak_info[0] <= $transcript_info[4] && $peak_info[0]>=$transcript_info[10] )
		{print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\t5'UTR"."\n";
		$type{"Peak_".$peak_info[0]."_5'UTR"}=1;
		}
		elsif ( $transcript_info[6] eq "+" && $peak_info[0] <= $transcript_info[5] && $peak_info[0]>= $transcript_info[11] )
                {
                                print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\t3'UTR"."\n";
			$type{"Peak_".$peak_info[0]."_3'UTR"}=1;
                }elsif ( $transcript_info[6]eq "-" && $peak_info[0] >= $transcript_info[5] && $peak_info[0] <= $transcript_info[11] )
                {print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\t3'UTR"."\n";
		$type{"Peak_".$peak_info[0]."_3'UTR"}=1;	
		}
		
		#Exon_intron
		else{
		my $ExCount=$transcript_info[7];
                my @ExStarts=split(",",$transcript_info[8]);
                my @ExEnds=split (",", $transcript_info[9]);

		for (my $i=0; $i<$ExCount;$i++)
                {
                if(($ExEnds[$i] >= $peak_info[0] && $peak_info[0]>= $ExStarts[$i]))
                {print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tExon"."\n";
		$type{"Peak_".$peak_info[0]."_Exon"}=1;
                last;}
                elsif( $i >= 1 && $ExEnds[($i-1)] <=$peak_info[0] && $ExStarts[$i]>=$peak_info[0])
                {
                print OUT "Peak_$peak_info[0]"."\t".$transcript_info[0]."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tIntron"."\n";
		$type{"Peak_".$peak_info[0]."_Intron"}=1;
                last;
                }
                }
	}
		




		}
		}

	}
}


my @states="";
my %count=();

foreach my $region (keys %type)
{
	@states=("5'UTR","3'UTR","Intron","Exon","Upstream:1000","Upstream:3000","Upstream:10000","Downstream:1000","Downstream:3000","Downstream:10000");
	for(my $i=0; $i< scalar(@states);$i++)
	{
		if ( $region =~/$states[$i]$/)
			{$count{$states[$i]}++;}
	}
		
}

open OUT_2, "+>","$ARGV[2].xls" or die "opening file: $ARGV[2].xls";

foreach my $state ( keys %count)
{
print OUT_2 $state."\t".$count{$state}."\n";
}



	

close PEAKS;
close GENES;
close OUT;
close OUT_2;




