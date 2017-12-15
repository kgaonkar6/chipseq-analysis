#!/usr/bin/perl
use strict;
use warnings;

if ( scalar(@ARGV) != 4)
{
print "./generate_report.pl wandy_BinCovg.txt SV.txt output_name output_dir";
exit 1;
}

my $BinCovg=$ARGV[0];
my $SV=$ARGV[1];
my $SV_name=`basename $SV .txt`;
chomp ($SV_name);
my $WORK_DIR=$ARGV[3];
my $output_name=$ARGV[2];
my $output=$WORK_DIR."/".$output_name.".conf";
open OUT,">","$output" or die "opening file: $output";
open IN ,"$BinCovg" or die "opening file: $BinCovg";
open IN_POS,">",$WORK_DIR."/".${output_name}."_pos_CNV.txt" or die "opening file ".$WORK_DIR."/".${output_name}."_pos_CNV.txt"; 
open IN_NEG,">",$WORK_DIR."/".${output_name}."_neg_CNV.txt" or die "opening file ".$WORK_DIR."/".${output_name}."_neg_CNV.txt";
open IN_NORM, ">" , $WORK_DIR."/".${output_name}."_norm_CNV.txt" or die "opening file".$WORK_DIR."/".${output_name}."_norm_CNV.txt";


my @valueCNV;
my $avg_norm_covg;
my %data_avg;

while (<IN>) {
	my $row=$_;
	chomp ($row);
	my @line=split("\t",$row);
	push @valueCNV,$line[3];
	$data_avg{"$line[0]\t$line[1]\t$line[2]"}="$line[3]";
	}

$avg_norm_covg=avg(@valueCNV);
print $avg_norm_covg;


foreach my $data (keys %data_avg) {
print IN_NORM "$data\t$avg_norm_covg\n";
our $value=$data_avg{$data}-$avg_norm_covg;
if ( $value >= 0 ){
print IN_POS "$data\t$value\n";}
if ( $value <0 ){
print IN_NEG "$data\t$value\n";}

}


close IN;



sub avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
         return $total / @_;
        }


print OUT "
karyotype=/data2/bsi/RandD/MAPRSeq/2.0/refs/tophat/karyotype.human.hg19.txt
show_ticks=yes
show_tick_lables=yes
<plots>
<plot>

show    = yes
type             = scatter
color            = green
file             = $WORK_DIR/${output_name}_pos_CNV.txt


r0 = 1.2r
r1 = 1.6r
min = 0
max = 4
glyph   = circle
glyph_size      = 10
stroke_color     = dgreen
stroke_thickness = 1
orientation=out
</plot>

<plot>

show    = yes
type             = scatter
color            = red
file             = $WORK_DIR/${output_name}_neg_CNV.txt


r0 = 1.0r
r1 = 1.2r
min = -4
max = 0
glyph   = circle
glyph_size      = 10
stroke_color     = dred
stroke_thickness = 1
orientation=out
</plot>



<plot>
show=yes
type=scatter
color =blue
file=$WORK_DIR/${output_name}_norm_CNV.txt
r0 = 1.2r
r1 = 1.2r
min = 0
max = 4
z=10
glyph   = circle
glyph_size      = 3
stroke_color     = blue
stroke_thickness = 3
orientation=out
</plot>

</plots>

<links>
z =0
color         = dblue
radius        = 0.6r
bezier_radius = 0.1r

<link SV>
show = yes
file = $SV
thickness     = 0.05
record_limit = 2500
</link>
</links>

<plot>
type=text
color=black
file=$WORK_DIR/${SV_name}_gene.txt
r0=0.6r
r1=0.7r
show_links=no
show_lable=yes
label_size=32p
label_font=condensed

padding=0p
rpadding=0p
</plot>



<ideogram>

<spacing>
default = 0.01r
break = 0.5r
</spacing>

radius    = 0.75r
thickness = 30p
fill=yes

stroke_color     = black
stroke_thickness = 2p

show_label       = yes
label_font       = condensed
label_radius     = 1r-75p
label_size       = 32p
label_parallel   = yes


</ideogram>
################################################################
## The remaining content is standard and required. It is imported 
## from default files in the Circos distribution.
##
## These should be present in every Circos configuration file and
## overridden as required. To see the content of these files, 
## look in etc/ in the Circos distribution.
#
<image>
## Included from Circos distribution.
<<include etc/image.conf>>
</image>
#
## RGB/HSV color definitions, color lists, location of fonts, fill patterns.
## Included from Circos distribution.
<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
</colors>
#
<fonts>
<<include etc/fonts.conf>>
</fonts>
#
## Debugging, I/O an dother system parameters
## Included from Circos distribution.
<<include etc/housekeeping.conf>>
chromosomes_units = 1000000
chromosomes_display_default = yes

";


close OUT;

#needs PERL5LIB=$PERL5LIB:/data2/bsi/reference/perl_workflow_ref/lib/perl5/x86_64-linux/auto:/usr/local/biotools/perl/5.10.0/bin/cpan:/usr/local/biotools/perl/5.10.0/bin/perl

`/projects/bsi/bictools/apps/visualizers/circos/circos-0.54/bin/circos -conf $WORK_DIR/$output_name.conf -png -outputfile $output_name.png -outputdir $WORK_DIR`;

