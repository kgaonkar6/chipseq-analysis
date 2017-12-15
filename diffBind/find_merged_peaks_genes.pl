#!/usr/bin/perl
use strict;
use warnings;
## Perl script takes in peak positions and finds genes with nearby TSS or TTS and reports their distance
## Jared Evans 12/19/11
## modified Krutika Gaonkar 12/18/15
use Data::Dumper;

if (scalar(@ARGV) != 6)	
{
	die ( "USAGE: find_nearby_genes.pl [peaks BED file] [reflat file] [output file] [up dist] [down dist] [merged bed] \n" );
}

open PEAKS, "$ARGV[0]" or die "opening file: $ARGV[0]";
open GENES, "$ARGV[1]" or die "opening file: $ARGV[1]";
open OUT, "+>","$ARGV[2]" or die "opening file: $ARGV[2]";

my $up_dist = $ARGV[3];
my $down_dist= $ARGV[4];
open MERGE,"$ARGV[5]" or die "opening file: $ARGV[0]";
my %peaks = ();
## chr_start_stop -> [peakcenter, FDR, chr, start, stop, peaklength]
my %transcripts = ();
## chr -> transcript_start_stop -> [transcript, gene, transcript start, transcript stop, TSS, TTS, strand,ExonCount,ExStarts,ExEnds,cdsStart,cdsEnd]
my %peaks_and_transcripts = ();
## peakchr_start_stop -> transcript_start_stop -> [transcript, gene, chr, start, stop, strand, TSS, distupstream_TSS, distdownstream_TSS, TTS, distupstream_TTS, distdownstream_TTS]

#my %type=();
## unique values

while (<MERGE>){
	my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
	$peaks_and_transcripts{"$line[0]_$line[1]_$line[2]"}=[$line[0],$line[1],$line[2],"NA","NA","NA","NA","NA","NA","NA","NA"];
}




open OUT, "+>","$ARGV[2]" or die "opening file: $ARGV[2]";
# load peaks into hash
while (<PEAKS>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $peakcenter = (($line[2]-$line[1])/2)+$line[1];
	if($peakcenter =~ m/\.5/){
		$peakcenter += 0.5;
	}
	$peaks{"$line[0]_$line[1]_$line[2]"} = [$peakcenter,$line[8],$line[0],$line[1],$line[2],($line[2]-$line[1]),$line[6],$line[7]];
	$peaks_and_transcripts{"$line[0]_$line[1]_$line[2]"}=[$line[0],$line[1],$line[2],"NA","NA","NA","NA","NA",$line[6],$line[7],$line[8]];
}
	#load all peaks once in output hash


# load transcripts into hash of hashes
while(<GENES>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $tss = 0;
	my $tts = 0;
	my $cdsStart=0;
	my $cdsEnd=0;
	my $strand = 1;
	if($line[3] eq "+"){
		$tss = $line[4];
		$tts = $line[5];
		$cdsStart= $line[6];
		$cdsEnd= $line[7];
	}else{
		$tss = $line[5];
		$tts = $line[4];
		$cdsStart=$line[7];
		$cdsEnd=$line[6];
	}
	$transcripts{$line[2]}{"$line[1]_$line[4]_$line[5]"} = [$line[1],$line[0],$line[4],$line[5],$tss,$tts,$line[3],$line[8],$line[9],$line[10],$cdsStart,$cdsEnd];
}

print OUT "Chr\tStart\tStop\tGene\tTSS\tTTS\tStrand\tINFO\tFC\tpvalue\tFDR\tDiffBind\n";

	
foreach my $peak_line (keys %peaks){
	#my $peak_line = "chr1_5011760_5014166"; #upstream -
	#my $peak_line = "chr1_15786175_15786773"; #upstream +
	#my $peak_line = "chr1_13355099_13355685"; #downstream -
	#my $peak_line = "chr1_34076808_34079752"; #downstream +
	my @peak_info = @{$peaks{$peak_line}};
	foreach my $transcript (keys %{$transcripts{$peak_info[2]}}){
		our @transcript_info = @{$transcripts{$peak_info[2]}{$transcript}};
		## TSS:
		if(($peak_info[0]-$transcript_info[4]) <= $down_dist && ($peak_info[0]-$transcript_info[4]) >=0 && $transcript_info[6] eq "+"){
#				print OUT "Peak_$peak_info[0]"."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tDownstream_TSS_$down_dist"."\t".$peak_info[6]."\t".$peak_info[7]."\t".$peak_info[1]."\n";
				$peaks_and_transcripts{"$peak_info[2]_$peak_info[3]_$peak_info[4]"}=[$peak_info[2],$peak_info[3],$peak_info[4],$transcript_info[1],$transcript_info[4],$transcript_info[5],$transcript_info[6],"Downstream_TSS_$down_dist",$peak_info[6],$peak_info[7],$peak_info[1]];
		}elsif(($transcript_info[4]-$peak_info[0]) <= $down_dist && ($transcript_info[4]-$peak_info[0]) >=0 && $transcript_info[6] eq "-"){
#				print OUT "Peak_$peak_info[0]"."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tDownstream_TSS_$down_dist"."\t".$peak_info[6]."\t".$peak_info[7]."\t".$peak_info[1]."\n";
				$peaks_and_transcripts{"$peak_info[2]_$peak_info[3]_$peak_info[4]"}=[$peak_info[2],$peak_info[3],$peak_info[4],$transcript_info[1],$transcript_info[4],$transcript_info[5],$transcript_info[6],"Upstream_TSS_$up_dist",$peak_info[6],$peak_info[7],$peak_info[1]];
		}
		elsif (($peak_info[0]-$transcript_info[4]) <= $up_dist && ($peak_info[0]-$transcript_info[4]) >=0 && $transcript_info[6] eq "-"){
  #                              print OUT "Peak_$peak_info[0]"."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tUpstream_TSS_$up_dist"."\t".$peak_info[6]."\t".$peak_info[7]."\t".$peak_info[1]."\n";
				$peaks_and_transcripts{"$peak_info[2]_$peak_info[3]_$peak_info[4]"}=[$peak_info[2],$peak_info[3],$peak_info[4],$transcript_info[1],$transcript_info[4],$transcript_info[5],$transcript_info[6],"Downstream_TSS_$down_dist",$peak_info[6],$peak_info[7],$peak_info[1]];
		}
		elsif(($transcript_info[4]-$peak_info[0]) <= $up_dist && ($transcript_info[4]-$peak_info[0]) >=0 && $transcript_info[6] eq "+"){
#                           print "Peak_$peak_info[0]"."\t".$peak_info[2]."\t".$transcript_info[1]."\t".$transcript_info[2]."\t".$transcript_info[3]."\t".$transcript_info[4]."\t".$transcript_info[5]."\t".$transcript_info[6]."\tUpstream_TSS_$up_dist"."\t".$peak_info[6]."\t".$peak_info[7]."\t".$peak_info[1]."\n";
				$peaks_and_transcripts{"$peak_info[2]_$peak_info[3]_$peak_info[4]"}=[$peak_info[2],$peak_info[3],$peak_info[4],$transcript_info[1],$transcript_info[4],$transcript_info[5],$transcript_info[6],"Upstream_TSS_$up_dist",$peak_info[6],$peak_info[7],$peak_info[1]];
		}

		}
		

	}

foreach my $peak_gene (sort keys %peaks_and_transcripts){
my @gene_info=@{$peaks_and_transcripts{$peak_gene}};
if ($gene_info[8] eq "NA")
{
print OUT $gene_info[0]."\t".$gene_info[1]."\t".$gene_info[2]."\t".$gene_info[3]."\t".$gene_info[4]."\t".$gene_info[5]."\t".$gene_info[6]."\t".$gene_info[7]."\t".$gene_info[8]."\t".$gene_info[9]."\t".$gene_info[10]."\t"."Not_significant"."\n";
}
elsif ($gene_info[8] > 0)
{
print OUT $gene_info[0]."\t".$gene_info[1]."\t".$gene_info[2]."\t".$gene_info[3]."\t".$gene_info[4]."\t".$gene_info[5]."\t".$gene_info[6]."\t".$gene_info[7]."\t".$gene_info[8]."\t".$gene_info[9]."\t".$gene_info[10]."\t"."Up"."\n";
}
elsif ($gene_info[8] < 0)
{
print OUT $gene_info[0]."\t".$gene_info[1]."\t".$gene_info[2]."\t".$gene_info[3]."\t".$gene_info[4]."\t".$gene_info[5]."\t".$gene_info[6]."\t".$gene_info[7]."\t".$gene_info[8]."\t".$gene_info[9]."\t".$gene_info[10]."\t"."Down"."\n";
}

}





	

close PEAKS;
close GENES;
close OUT;





