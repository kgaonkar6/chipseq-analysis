#!/usr/bin/perl


use strict;
use warnings;

use File::Basename;

$|=1;

# Expected input:
# -i path/projectName_sample1.segment.vcf, path/projectName_sample2.segment.vcf, ...




# --------------------------------------------------------------------------------
# declare global variables/DS
# --------------------------------------------------------------------------------

my ($vcfPathString, $outputFile, $displayUsage, @tp);
my $numFormatFields=0;
my $SCRIPTNAME="mergeSegmentVCFs.pl";




# --------------------------------------------------------------------------------
# process args
# --------------------------------------------------------------------------------

if ($#ARGV < 1) { warn "Insufficient args specified\n"; &usage; exit 1; }

$vcfPathString=$ARGV[0];
$outputFile=$ARGV[1];

if (!$vcfPathString     		){ warn "\nMust specify comma-separated list of .segment.vcf files to merge\n"; &usage; exit 1; }
if (!$outputFile       		){ warn "\nMust specify path to output merged VCF\n"; &usage; exit 1; }

my @vcfPaths=split(',',$vcfPathString);




# --------------------------------------------------------------------------------
# Write merged VCF to $outputFile
# --------------------------------------------------------------------------------


# get sample names

 my @sampleNames=();
 foreach my $vcfPath (@vcfPaths){
	my $vcfFilename=basename($vcfPath);
	@tp=split('_',$vcfFilename);
	my $prefix=$tp[0].'_';
	$vcfFilename =~ s/$prefix//;
	$vcfFilename =~ s/\.segment\.vcf//;
	push @sampleNames, $vcfFilename;
 }

 
# write all header fields before ##SAMPLE
 
open(OUT,'>',$outputFile) or die $!;

open(IN, '<', $vcfPaths[0]) or die $!;
while(<IN>){

		if($_ =~ /^##SAMPLE/){
			last;
		} else{
			print OUT $_;
		}

}
close(IN);


# print ##SAMPLE line from each vcf file

foreach my $i ( 0..(scalar(@vcfPaths)-1) ){

	my $vcfPath = $vcfPaths[$i];
	open(IN, '<', $vcfPath) or die $!;
	while(<IN>){
		if($_ =~ /^##SAMPLE/){
			print OUT $_;
			last;
		}	
	}

}
close(IN);


# write all header fields after ##SAMPLE

open(IN, '<', $vcfPaths[0]) or die $!;
my $write=0;
while(<IN>){
		if($write && $_ =~ /^##/){
			print OUT $_;
		} else{
			if($_ =~ /^##SAMPLE/){
				$write=1;
			} 
		}
}
close(IN);


# write rows

foreach my $i ( 0..(scalar(@vcfPaths)-1) ){

	my $vcfPath=$vcfPaths[$i];
	my $currentSampleName=basename($vcfPath);
	@tp=split('_', $currentSampleName);
	my $prefix=$tp[0].'_';
	$currentSampleName =~ s/$prefix//;
	$currentSampleName =~ s/\.segment\.vcf//;
	
	open(IN, '<', $vcfPath) or die $!;
	while(<IN>){
		chomp;
		my $line = $_;
		if($i == 0 && $line =~ /^##/){
			if($line =~ /^##FORMAT/){
				$numFormatFields = $numFormatFields + 1;
			}
		} elsif( $i == 0 && $line =~/^#CHROM/){
			print OUT join("\t",'#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', @sampleNames), "\n";
		} elsif( $line =~/^#/){
			next;
		}
		else{
			my @rowVals = split("\t",$line);
			foreach my $j ( 0..scalar(@rowVals)-2){
				print OUT $rowVals[$j],"\t";
			}
			foreach my $sampleIndex ( 0..$#sampleNames ){
				if($sampleNames[$sampleIndex] eq $currentSampleName){
					print OUT $rowVals[$#rowVals],"\t";
				} else{
					print OUT '.';
					print OUT ':.' x ($numFormatFields-1);
					print OUT "\t";
				}
			}
			print OUT "\n";
		}
	}	
	close(IN);
}

close(OUT);

print "\nMerged VCF written to $outputFile.\n\n";




# --------------------------------------------------------------------------------
# FUNCTIONS
# --------------------------------------------------------------------------------

sub usage { print "\n  Usage: $SCRIPTNAME <vcf1,vcf2,..> <outFile>\n\n";
            #print "\n\t\t -i vcf1,vcf2,..";
			#print "\n\t\t -o outfile..";
			#print "\n\t\t -h display usage\n\n"   ;
}




