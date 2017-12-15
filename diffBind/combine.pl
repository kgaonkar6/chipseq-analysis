#!/usr/bin/perl
use strict;
use warnings;
## Perl script takes in peak positions and finds genes with nearby TSS or TTS and reports their distance
## Jared Evans 12/19/11
## modified Krutika Gaonkar 12/18/15
use Data::Dumper;



if (scalar(@ARGV) != 5)	
{
	die ( "USAGE: find_nearby_genes.pl [peaks BED file] [reflat file] [output file (Mark_group.*xls] [up dist] [down dist] \n" );
}

open PEAKS, "$ARGV[0]" or die "opening file: $ARGV[0]";
open GENES, "$ARGV[1]" or die "opening file: $ARGV[1]";
open OUT, "+>","$ARGV[2]" or die "opening file: $ARGV[2]";

my $up_dist = $ARGV[3];
my $down_dist= $ARGV[4];
my %peaks = ();
## chr_start_stop -> [ peakcenter, chr, start, stop, FC, pvalue, FDR]
my %genes = ();
## chr_start_stop -> [ chr, cd_start, cd_stop, gene, gene_details, strand ]
my %peaks_and_genes = ();
## chr_start_stop_gene ->[ gene_chr, cds_start, cds_stop, gene, gene_details, strand, INFO, FC, pvalue, FDR, DB_INFO]

#my %type=();
## unique values
my $group=`basename $ARGV[2]|cut -d "." -f1`;
chomp($group);



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
	$peaks{"$line[0]_$line[1]_$line[2]"} = [$peakcenter,$line[0],$line[1],$line[2],$line[6],$line[7],$line[8]];
}


while(<GENES>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $cd_start=0;
	my $cd_stop=0;
	if($line[5] eq "+"){
		$cd_start = $line[1];
		$cd_stop = $line[2];
	}else{
		$cd_start = $line[2];
		$cd_stop = $line[1];
	}
	$genes{"$line[0]_$line[1]_$line[2]_$line[3]"} = [$line[0],$cd_start,$cd_stop,$line[3],$line[4],$line[5]];
	$peaks_and_genes{"$line[0]_$line[1]_$line[2]_$line[3]"}=[$line[0],$cd_start,$cd_stop,$line[3],$line[4],$line[5],"NA",0,0,"NA",0,0,0,"NA"];
}



print OUT "Chr\tStart\tStop\tGene\tDetails\tStrand\tpeak_Chr\tpeak_Start\tpeak_Stop\tINFO\t${group}_FC\t${group}_pvalue\t${group}_FDR\t${group}_DB\n";

	
foreach my $gene_region ( keys %peaks_and_genes){
	our @gene_info = @{$genes{$gene_region}};
	our $DB="NA";
	foreach my $peak_region ( keys %peaks){
		our @peak_info = @{$peaks{$peak_region}};
		if ($gene_info[0] eq $peak_info[0]){
		if ( $peak_info[4] == 0){$DB="NA";}
		elsif ($peak_info[4] >0){$DB="Up";}
		elsif ($peak_info[4] <0){$DB="Down";}
		
		## TSS:
		if(($peak_info[0]-$gene_info[1]) <= $down_dist && ($peak_info[0]-$gene_info[1]) >=0 && $gene_info[5] eq "+"){
				$peaks_and_genes{"$gene_info[0]_$gene_info[1]_$gene_info[2]_$gene_info[3]"}=[$gene_info[0],$gene_info[1],$gene_info[2],$gene_info[3],$gene_info[4],$gene_info[5],$peak_info[1],$peak_info[2],$peak_info[3],"Downstream_$down_dist",$peak_info[4],$peak_info[5],$peak_info[6],$DB];
		}
		elsif(($peak_info[0]-$gene_info[1]) <= $up_dist && ($peak_info[0]-$gene_info[1]) >=0 && $gene_info[5] eq "-"){
                                $peaks_and_genes{"$gene_info[0]_$gene_info[1]_$gene_info[2]_$gene_info[3]"}=[$gene_info[0],$gene_info[1],$gene_info[2],$gene_info[3],$gene_info[4],$gene_info[5],$peak_info[1],$peak_info[2],$peak_info[3],"Upstream_$down_dist",$peak_info[4],$peak_info[5],$peak_info[6],$DB];

		}
		elsif(($gene_info[1]-$peak_info[0]) <= $up_dist && ($gene_info[1]-$peak_info[0]) >=0 && $gene_info[5] eq "+"){
                               $peaks_and_genes{"$gene_info[0]_$gene_info[1]_$gene_info[2]_$gene_info[3]"}=[$gene_info[0],$gene_info[1],$gene_info[2],$gene_info[3],$gene_info[4],$gene_info[5],$peak_info[1],$peak_info[2],$peak_info[3],"Upstream_$up_dist",$peak_info[4],$peak_info[5],$peak_info[6],$DB];
		}
		elsif(($gene_info[1]-$peak_info[0]) <= $down_dist && ($gene_info[1]-$peak_info[0]) >=0 && $gene_info[5] eq "-"){
                               $peaks_and_genes{"$gene_info[0]_$gene_info[1]_$gene_info[2]_$gene_info[3]"}=[$gene_info[0],$gene_info[1],$gene_info[2],$gene_info[3],$gene_info[4],$gene_info[5],$peak_info[1],$peak_info[2],$peak_info[3],"Downstream_$down_dist",$peak_info[4],$peak_info[5],$peak_info[6],$DB];
		}	
	}
}	
	
}

foreach my $peak_gene ( keys %peaks_and_genes){
my @peak_gene_info=@{$peaks_and_genes{$peak_gene}};
print OUT $peak_gene_info[0]."\t".$peak_gene_info[1]."\t".$peak_gene_info[2]."\t".$peak_gene_info[3]."\t".$peak_gene_info[4]."\t".$peak_gene_info[5]."\t".$peak_gene_info[6]."\t".$peak_gene_info[7]."\t".$peak_gene_info[8]."\t".$peak_gene_info[9]."\t".$peak_gene_info[10]."\t".$peak_gene_info[11]."\t".$peak_gene_info[12]."\t".$peak_gene_info[13]."\n";


}





	

close PEAKS;
close GENES;
close OUT;





