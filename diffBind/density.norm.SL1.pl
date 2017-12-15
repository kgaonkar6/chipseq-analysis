#!/usr/bin/perl
use strict;

use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                                                  'runinfo|r=s',
                                                
                                                  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){

    print "Usage:\nruninfo|r=s\nmergedPeaks|mrgpk=s\nhelp|h";
    pod2usage( {-exitval => 1, -verbose => 2, -output => \*STDERR} );
            }
my $run_info=$options{runinfo};
my $WORK_DIR=`grep -w "WORK_DIR" $run_info|cut -d "=" -f2`;
my $BEDTOOLS=`grep -w "BEDTOOLS" $run_info|cut -d "=" -f2`;
my $WORKDIR=`grep -w "WORKDIR" $run_info|cut -d "=" -f2`;

opendir (I_DIR, $WORKDIR) or die $!;


while ( my $file = readdir (I_DIR)){
        $file =~ s/\-/\_/g;
        if ($file =~ /[Ii]nput.+\S+density.bed$/){
        @sample_base= split(/_[Ii]nput/, $file);
	my $treat_IP=`ls *$sample_base[0]*density.bed`;
	chomp ($treat_IP);
	my $treat_input=$file;
	chomp ($treat_input);
	
	}
}
close OUT;

