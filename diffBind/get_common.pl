#!/usr/bin/perl
use strict;
use warnings;
## Perl script takes DB bed files to get occupancy for K4me3/K27me3/K9me3

use Data::Dumper;



if (scalar(@ARGV) != 2)	
{
	die ( "USAGE: combine.pl [DB peaks BED file] [output file]  \n" );
}

open GENES, "$ARGV[0]" or die "opening file: $ARGV[0]";
open OUT, "+>","$ARGV[1]" or die "opening file: $ARGV[1]";

my %genes =();
my %peaks_and_genes = ();
## chr_start_stop_gene ->[ gene_chr, cds_start, cds_stop, gene, gene_details, strand, INFO, FC, pvalue, FDR]




open OUT, "+>","$ARGV[1]" or die "opening file: $ARGV[1]";



while(<GENES>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	#chr1	11874	14408	DDX11L1	DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1	+ 
	$peaks_and_genes{"$line[3]"} =[$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[15],$line[16],$line[17]];
	$genes{"$line[3]"} = [$line[0],$line[1],$line[2],$line[3],$line[4],$line[5],"NA","NA","NA","NA","NA","NA"];
	
}



print OUT "Chr\tStart\tStop\tGene\tDetails\tStrand\tK4me3\tK27me3\tK9me3\tK4me3+27me3\tK4me3+K9me3\tK27me3+K9me3\n";

	

foreach my $gene_region ( keys %genes){
	my @gene_info = @{$genes{$gene_region}};
	my @peak_gene_info=@{$peaks_and_genes{"$gene_region"}};
	my $K4me3=$gene_info[6];
	my $K27me3=$gene_info[7];
 	my $K9me3=$gene_info[8];
	my $K4me3_K27me3=$gene_info[9];
	my $K4me3_K9me3=$gene_info[10];
	my $K27me3_K9me3=$gene_info[11];
	if ($peak_gene_info[0]!=0||$peak_gene_info[1]!=0||$peak_gene_info[2]!=0	){
	$K4me3="Y";
	}
	if ($peak_gene_info[3]!=0||$peak_gene_info[4]!=0||$peak_gene_info[5]!=0 ){
      	$K27me3="Y";
	if ($K4me3 eq "Y"){$K4me3_K27me3="Y";}
        }
	if ($peak_gene_info[6]!=0||$peak_gene_info[7]!=0||$peak_gene_info[8]!=0||$peak_gene_info[9]!=0||$peak_gene_info[10]!=0||$peak_gene_info[11]!=0){
        $K9me3="Y";
	if ($K4me3 eq "Y"){$K4me3_K9me3="Y";}
	if ($K27me3 eq "Y"){$K27me3_K9me3="Y";}
	}
	$genes{"$gene_info[3]"} = [$gene_info[0],$gene_info[1],$gene_info[2],$gene_info[3],$gene_info[4],$gene_info[5],$K4me3,$K27me3,$K9me3,$K4me3_K27me3,$K4me3_K9me3,$K27me3_K9me3];
	
}

foreach my $gene (sort keys %genes){
my @gene_info=@{$genes{$gene}};
print OUT $gene_info[0]."\t".$gene_info[1]."\t".$gene_info[2]."\t".$gene_info[3]."\t".$gene_info[4]."\t".$gene_info[5]."\t".$gene_info[6]."\t".$gene_info[7]."\t".$gene_info[8]."\t".$gene_info[9]."\t".$gene_info[10]."\t".$gene_info[11]."\n";

}





	

close GENES;
close OUT;





