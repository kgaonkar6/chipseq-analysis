#!/usr/bin/perl

# This wrapper script calls PatternCNV to perform tertiary CNV analyses on secondary analysis bams against a list of reference samples. 


use 5.010001;

use File::Basename;
use File::Copy;
use File::Spec;
use Getopt::Long;

use strict;
use warnings;

$|=1;




# --------------------------------------------------------------------------------
# declare global variables/DS
# --------------------------------------------------------------------------------

my $SCRIPTNAME="run_cnvAnalyses.pl";

my $JOBLIMIT  = 0; 
my $RUNTIMELIMIT = 10000000;

my $sub_w = `whoami`; chomp ( $sub_w );
my $num_jobs_running = 1;

my ($indir, $outdir, $project, $email, $panelconfig, $mainconfig, $sampleMetaData, $queue, $stageExpr, $debug, $regenWigs, $displayUsage);

my %sample_gender = ();
my $genderKnown = 0;

my %cfg_hash = ();
my @SAMPLEIDS = ();
my @GENEIDS = ();
my @jobids=();
my $stdout='';
my $jobIDString='';
my $holdJID=' ';
my $exec='';

my $time = localtime;

my @chars = ();
my $i = 0; # change to something more meaningful

my $errors=0;

my @tp=();




# --------------------------------------------------------------------------------
# process args
# --------------------------------------------------------------------------------

if ($#ARGV < 4) { warn "Insufficient args specified\n"; &usage; exit 1; }

my @Options =
(
    "id=s"  => \$indir,
    "od=s"  => \$outdir,
    "pj=s"  => \$project,
    "mc=s"  => \$mainconfig,
    "pc=s"  => \$panelconfig,
    "sd=s"  => \$sampleMetaData,
    "rw=s"  => \$regenWigs,
    "q=s"   => \$queue,
    "e=s"   => \$email,
    "s=s"   => \$stageExpr,
    "dg=s" => \$debug,
	"h" => \$displayUsage
);


if (!&GetOptions(@Options) ){ &usage; exit 1;}
if ( defined $displayUsage ){ &usage; exit 1; }
if (!$indir       		){ warn "\nMust specify an input directory\n"; &usage; exit 1; } $indir =~ s/\/$//; $indir =`cd "$indir"; pwd`; chomp $indir;
if (!$outdir       		){ warn "\nMust specify an output directory\n"; &usage; exit 1; } $outdir =~ s/\/$//; $outdir = "${outdir}/cnv_usingTrainingSet"; $outdir =`mkdir -p "$outdir"; cd "$outdir"; pwd`; chomp $outdir;
if (!$project     		){ warn "\nMust specify project name  \n"; &usage; exit 1; }
if (!$mainconfig  ){ warn "\nMust specify main config file\n"; &usage; exit 1; }
if (!$panelconfig       ){ warn "\nMust specify a panel config file \n"; &usage; exit 1; }
if (!$sampleMetaData       ){ warn "\nNote: Sample metadata file was not specified.\n"; }
if (!$debug       ){ $debug = 0; }
if (! defined $stageExpr ){ print "\n* No stage(s) specified. Using -s a"; $stageExpr = 'a'; }
my @stages = split /[\s,\,]+/, $stageExpr; 
@stages = map { lc } @stages;    
if ( ! ( ('a' ~~ @stages) || ('1' ~~ @stages) || ('2' ~~ @stages) || ('3' ~~ @stages) || ('4' ~~ @stages) ) ) { warn "\nInvalid stage number(s)\n"; &usage; exit 1; }
if ( 'a' ~~ @stages || (! $stageExpr) ) { @stages = qw( 1 2 3 4 ); print "\n* Note: Running all stages"; };
if (!$queue       ){ $queue       = "sandbox.q"; print "\n* No queue specified. Using -q $queue";}
if (!$email       ){ $email = `finger \$\(whoami\) | cut -f2 -d \'\;\' | head -n 1`; $email =~ s/\s+//g; print "\n* No email specified. Using -e $email"; }
if (! defined $regenWigs       ){  $regenWigs = 0; } else { print "\n* Wig files for the training samples will be regenerated"; $regenWigs = 1; }




# --------------------------------------------------------------------------------
# create log file
# --------------------------------------------------------------------------------

`mkdir -p ${outdir}/logs`;
`mkdir -p ${outdir}/logs/jobs`;
open ( LOG, "> ${outdir}/logs/main.log" ) or die "\nERROR: Unable to create log file.\n\n";

		


# --------------------------------------------------------------------------------
# check input
# --------------------------------------------------------------------------------

$time = localtime;
#print "\n\n[$time] Checking input...";
&printToLogNStdout("\n\n[$time] Checking input...");
$errors=0;

                                                 
if(! -d $indir){&printToLogNStdout( "ERROR: Invalid input directory $indir"); $errors += 1; }

opendir(my $DIR, $indir) or die $!;
#my @bamFiles = glob "*.bam";
my @bamFiles =`find $indir -name '*.bam'`;

print "\n Bamfiles: ", @bamFiles;

#while ( my $entry = readdir $DIR ) {
#    next unless -d $indir . '/' . $entry;
#    next if $entry eq '.' or $entry eq '..' or $entry eq 'lane1';
#    push(@SAMPLEIDS, $entry);
#}


for(0..$#bamFiles){
    my $sampleID = $bamFiles[$_];
    my $ext = '.bam';
$sampleID = basename($sampleID);
    $sampleID =~ s/$ext//g;
chomp($sampleID);
    push(@SAMPLEIDS, $sampleID);
}
closedir $DIR;

print " SampleIDs: ", join(',', @SAMPLEIDS);

if(@SAMPLEIDS < 2){	&printToLogNStdout( "ERROR: Insufficient samples (Min. required = 2)"); $errors += 1; }

#foreach my $sampleID (@SAMPLEIDS){
#	if(! (-f "${indir}/${sampleID}/${sampleID}.bam" ) ){
#		&printToLogNStdout( "\nERROR: Cannot find ${indir}/${sampleID}/${sampleID}.bam"); 
#		$errors += 1;
#	}
#	elsif( ! (-s "${indir}/${sampleID}/${sampleID}.bam" ) ){
#		&printToLogNStdout( "\nERROR: ${indir}/${sampleID}/${sampleID}.bam has zero size"); 
#		$errors += 1;
#	}
#}

if(-d $outdir){ chdir $outdir or &printToLogNStdout( "ERROR: Unable to enter $outdir") and $errors += 1; }

if( ! (-f -s $mainconfig )){ &printToLogNStdout( "ERROR: Invalid mainconfig file specified: $mainconfig"); $errors += 1; }

if( ! (-f -s $panelconfig )){ &printToLogNStdout( "ERROR: Invalid panelconfig file specified: $panelconfig"); $errors += 1; }

if( $sampleMetaData && (! (-f -s $sampleMetaData )) ){ &printToLogNStdout( "ERROR: Invalid sample metadata file specified: $sampleMetaData"); $errors += 1; }

&quit(2) if($errors);




# --------------------------------------------------------------------------------
# read in panelconfig
# --------------------------------------------------------------------------------
# example .cfg file: CAPN CANCP CLC kt /dlmp/sandbox/cgslIS/bin/scripts_clc6 /dlmp/sandbox/cgslIS/proj/BenKipp/bin/CANCP_NGS05.cfg 

`mkdir -p ${outdir}/configs`;

$time = localtime;
&printToLogNStdout( "\n\n[$time] Loading panel config ${panelconfig}...");
$errors=0;

my @reqdPanelConfigVars=qw( TRAININGSET_BAMS CANDIDATE_TRAININGSET_BAMS GENOME_REF GENOME_REF_NAME GENOME_SIZE EXON_BED CAPTUREKIT_BED COVG_RES BIN_SIZE MIN_MAPPING_QUAL EXTENSION_BUFFER MIN_NUM_SAMPLES MAX_MAD_LG2RATIO_DIFF MAX_STDEV_LG2RATIO MIN_CORR MIN_COVG MAX_NUM_DELS MAX_NUM_DUPS MIN_INPUT_CONC MIN_FRAGMENT_SIZE CUTOFF_LG2R_DUPS CUTOFF_LG2R_DELS CUTOFF_ZSCORE_DUPS CUTOFF_ZSCORE_DELS MIN_BINS QSUB_EXONKEY_MEMORY QSUB_BAM2WIG_MEMORY QSUB_CALLCNVS_MEMORY GAINS_HIGH_TRANCHE_MIN_PVAL GAINS_HIGH_TRANCHE_MAX_PVAL GAINS_MEDIUM_TRANCHE_MIN_PVAL GAINS_MEDIUM_TRANCHE_MAX_PVAL GAINS_LOW_TRANCHE_MIN_PVAL GAINS_LOW_TRANCHE_MAX_PVAL LOSSES_HIGH_TRANCHE_MIN_PVAL LOSSES_HIGH_TRANCHE_MAX_PVAL LOSSES_MEDIUM_TRANCHE_MIN_PVAL LOSSES_MEDIUM_TRANCHE_MAX_PVAL LOSSES_LOW_TRANCHE_MIN_PVAL LOSSES_LOW_TRANCHE_MAX_PVAL );


open ( IN, "< $panelconfig " ) or die "\n$!\n\n";
while ( <IN> ) { 
  if ( !/^\#/ ) {
    chomp; 
    my @tp = split ( /=/, $_, 2 );
    if(scalar(@tp) == 2){
    	if ( $tp[0] && $tp[1] ) { $cfg_hash{$tp[0]} = $tp[1]; }
    	elsif ( ! $tp[1] ) { &printToLogNStdout( "\n\tWarning: $tp[0] is not set"); }
    }
  }
}
close ( IN );

foreach my $reqdPanelConfigVar (@reqdPanelConfigVars){
	if(! $cfg_hash{$reqdPanelConfigVar}){
		&printToLogNStdout( "\n  Error: $reqdPanelConfigVar is a required parameter that has not been set. ");
		$errors += 1;
	}
}

&quit(2) if($errors);

copy($panelconfig, "${outdir}/configs") or die "$!";




# --------------------------------------------------------------------------------------------------------------------------------------
# read in mainconfig   Note: In the case of duplicate variables, mainconfig's vars will overwrite panelconfig's
# --------------------------------------------------------------------------------------------------------------------------------------
$time = localtime;
&printToLogNStdout( "\n\n[$time] Loading main config ${mainconfig}...");
$errors=0;

my @reqdMainConfigVars=qw( PATTERNCNV SAMTOOLS BEDTOOLS R RLIB PIPELINE QSTAT QSUB NJOBS_LIMIT );

open ( IN, "< $mainconfig " ) or die "\n$!\n\n";
while ( <IN> ) { 
  if ( !/^\#/ ) {
    chomp; 
    my @tp = split ( /=/, $_, 2 );

    if ( ( defined $tp[0] ) && ( defined $tp[1] ) ) { $cfg_hash{$tp[0]} = $tp[1]; }
    #elsif( ! $tp[0] ~~ @reqdMainConfigVars ) { &printToLogNStdout( "\n  Warning: $tp[0] is not a required parameter and will not be read."); }
    
  }
}
close ( IN );

foreach my $reqdMainConfigVar (@reqdMainConfigVars){
	if(! $cfg_hash{$reqdMainConfigVar}){
		&printToLogNStdout( "\n  Error: $reqdMainConfigVar is a required parameter that has not been set. ");
		$errors += 1;
	}
}

&quit(2) if($errors);
copy($mainconfig, "${outdir}/configs/main.cfg") or die "$!";




# --------------------------------------------------------------------------------
# check pipeline installation
# --------------------------------------------------------------------------------
$time = localtime;
&printToLogNStdout( "\n\n[$time] Checking pipeline installation...");
$errors=0;

my @reqdScripts=qw( create_cnv_txt.R gen.pVal.matrix.R gen.QC.metrics.wrtTrainingSet.R gen.seg.pVal.matrix.R gen.VCF_exon.R gen.VCF_segment.R );

foreach my $reqdScript (@reqdScripts){
	if( ! -f "$cfg_hash{PIPELINE}/scripts/${reqdScript}" ){
		&printToLogNStdout( "\n Error: Required script $reqdScript is missing in $cfg_hash{PIPELINE}/scripts");
		$errors += 1;
	}
}

&quit(2) if ($errors);




# --------------------------------------------------------------------------------
# QC sample metadata file
# --------------------------------------------------------------------------------
if ($sampleMetaData){

	$time = localtime;
	&printToLogNStdout( "\n\n[$time] Checking sample info file...");
	$errors=0;

	`mkdir -p ${outdir}/qc`;

	my @colNames = ();
	open(SAMPLE_METADATA, "< $sampleMetaData") or die "$!";
	open(SAMPLE_METADATA_QC, ">${outdir}/qc/metadata.qc.csv") or die "$!";

	while(<SAMPLE_METADATA>){

	 chomp;

	 if($. == 1){

	 	@colNames = map(uc,split(/\t/, $_));

	 	print "\n", join(',', @colNames), "\n";

	 	print SAMPLE_METADATA_QC join(',', @colNames), "\n";

	 }else{
	
	 	my @rowValues = split(/\t/, $_);

		print join(',', @rowValues); print "\n";

	 	my $sampleName = uc($rowValues[0]);

	 	for(my $colIndex = 0; $colIndex <= $#colNames; $colIndex++){
	

			#if($sampleName ~~ @SAMPLEIDS){

				$sample_gender{$sampleName} = uc($rowValues[$colIndex]) if ( lc($colNames[$colIndex]) eq "gender" || lc($colNames[$colIndex]) eq "sex");

				print SAMPLE_METADATA_QC "$rowValues[$colIndex]";
				print SAMPLE_METADATA_QC " (QC Fail)" if( defined $cfg_hash{"MIN_".uc($colNames[$colIndex])} && $rowValues[$colIndex] < $cfg_hash{"MIN_".uc($colNames[$colIndex])} && ($colIndex != 0 && ( lc($colNames[$colIndex]) ne "gender" || lc($colNames[$colIndex]) ne "sex" ) ));
				print SAMPLE_METADATA_QC " (QC Fail)" if( defined $cfg_hash{"MAX_".uc($colNames[$colIndex])} && $rowValues[$colIndex] > $cfg_hash{"MAX_".uc($colNames[$colIndex])} && ($colIndex != 0 && ( lc($colNames[$colIndex]) ne "gender" || lc($colNames[$colIndex]) ne "sex" ) ));
				print SAMPLE_METADATA_QC "," if $colIndex != $#colNames; 

			#}

	
   		}

   		print SAMPLE_METADATA_QC "\n";

		}
	}

	
if( 'SEX' ~~ @colNames || 'GENDER' ~~ @colNames ){
	$genderKnown = 1;
	print("\nGender column detected - PCNV will be run in gender-aware mode.");
}


# check if gender info is avail for target samples

 if($genderKnown){

  foreach my $sampleID (@SAMPLEIDS){

	if(! $sample_gender{$sampleID}){

		print "\n  Error: Gender information has not been specified for sample $sampleID";
		$errors = $errors + 1; 
	}	

  }

 }

close(SAMPLE_METADATA_QC);
close(SAMPLE_METADATA);
 
&quit(2) if($errors);


} else{

	$sampleMetaData="NA";

}




# --------------------------------------------------------------------------------
# generate exon key file
# --------------------------------------------------------------------------------
if($regenWigs){
$time = localtime;
&printToLogNStdout( "\n\n[$time] Generating exon key file...");
$errors=0;

`mkdir -p ${outdir}/pcnv/configs`;

# create config file to generate exon key

 open(EXON_KEY_CFG, "> ${outdir}/pcnv/configs/exon_key.sh.cfg") or die "\n$!\n\n";
 print EXON_KEY_CFG "# Paths to local installs of SAMtools and BEDtools";
 print EXON_KEY_CFG "\nPATTERNCNV=$cfg_hash{PATTERNCNV}";
 print EXON_KEY_CFG "\nSAMTOOLS=$cfg_hash{SAMTOOLS}";
 print EXON_KEY_CFG "\nBEDTOOLS=$cfg_hash{BEDTOOLS}";
 print EXON_KEY_CFG "\n\nGENOME_SIZE=$cfg_hash{GENOME_SIZE}";
 print EXON_KEY_CFG "\nGENOME_REF=$cfg_hash{GENOME_REF}";
 close(EXON_KEY_CFG);

# generate exon key file

 #print "$cfg_hash{PATTERNCNV}/bam2wig/exon_key.sh -e $cfg_hash{EXON_BED} -c $cfg_hash{CAPTUREKIT_BED} -o ${outdir}/pcnv/configs/exon_key.txt -b $cfg_hash{COVG_RES} -t ${outdir}/pcnv/configs/exon_key.sh.cfg -x $cfg_hash{EXTENSION_BUFFER} -s $cfg_hash{BIN_SIZE}";
 $exec="$cfg_hash{PATTERNCNV}/bam2wig/exon_key.sh -e $cfg_hash{EXON_BED} -c $cfg_hash{CAPTUREKIT_BED} -o ${outdir}/pcnv/configs/exon_key.txt -b $cfg_hash{COVG_RES} -t ${outdir}/pcnv/configs/exon_key.sh.cfg -x $cfg_hash{EXTENSION_BUFFER} -s $cfg_hash{BIN_SIZE}";
 &printToLogNStdout("\n\t$exec");
 #`$cfg_hash{PATTERNCNV}/bam2wig/exon_key.sh -e $cfg_hash{EXON_BED} -c $cfg_hash{CAPTUREKIT_BED} -o ${outdir}/pcnv/configs/exon_key.txt -b $cfg_hash{COVG_RES} -t ${outdir}/pcnv/configs/exon_key.sh.cfg -x $cfg_hash{EXTENSION_BUFFER} -s $cfg_hash{BIN_SIZE}`;
 $stdout=`echo "$exec" | $cfg_hash{QSUB} -N GEN.EXONKEY -m abe -M ${email} -V -wd ${outdir}/pcnv/logs -q $queue $cfg_hash{QSUB_EXONKEY_MEMORY}`;
 &printToLogNStdout("\n\t$stdout");
 @tp = split(' ', $stdout);
 $holdJID = $tp[2];

#&quit(2) if (! -f -s "${outdir}/pcnv/configs/exon_key.txt");
}



# --------------------------------------------------------------------------------
# check training set BAMs and convert to wig if absent 
# --------------------------------------------------------------------------------
$time = localtime;
&printToLogNStdout( "\n\n[$time] Checking training BAMs...");
$errors=0;
@jobids = ();

`mkdir -p ${outdir}/pcnv/wigs`;
`mkdir -p ${outdir}/pcnv/logs`;


# Remove breadcrumb files if present
 unlink "${outdir}/pcnv/sortTrainingBam.complete";


open ( IN, "< $cfg_hash{TRAININGSET_BAMS}" ) or die "\n$!\n\n";
while(<IN>){

	chomp;
	my $sortedBamPath = $_;
	my $bamDir = dirname($sortedBamPath);
	#my $bam = basename($sortedBamPath);
	my $sortedBam = basename($sortedBamPath);
	chdir($bamDir) or print "Cannot change to $sortedBam dir";
	#my $sortedBam = $bam;
	#$sortedBam =~ s/\.bam/\.sorted\.bam/;
	my $wigFile = "${sortedBam}.coverage.wig.gz";
	#print "\n\n${wigFile}\n\n";

	if( (! -s -f $wigFile) || $regenWigs ){
		
		if(! -s -f $sortedBam ){

			##todo

			if(-l $sortedBam){
				
			 	if (lstat $sortedBam and not stat $sortedBam){
					print "\n Error: Symlink to sorted training BAM $sortedBam is broken"; 
					$errors = $errors + 1;
				} 

			} else{

				# sort bam

				print "\n\tError: Training BAM $sortedBamPath is non-existent/empty."; 
				$errors = $errors + 1;

			}

					
		} else{

			&printToLogNStdout("\n\t* Launching wig-gen job for $sortedBam...");
			my $stdout = `echo \"$cfg_hash{"PATTERNCNV"}/bam2wig/bam2wig.sh -i $sortedBamPath -o $bamDir -b $cfg_hash{"COVG_RES"} -m $cfg_hash{"MIN_MAPPING_QUAL"} -t ${outdir}/pcnv/configs/exon_key.sh.cfg -e ${outdir}/pcnv/configs/exon_key.txt; cp ${bamDir}/${sortedBam}.coverage.wig.gz ${outdir}/pcnv/wigs\" \| qsub -N TRAINING.BAM2WIG_${sortedBam} -V -wd ${outdir}/pcnv/logs -q $queue $cfg_hash{"QSUB_BAM2WIG_MEMORY"}`;
			&printToLogNStdout("\n\t$stdout");
			my @tp = split(' ', $stdout);
			push(@jobids, $tp[2]);

		}


	} else{

		copy("${sortedBam}.coverage.wig.gz", "${outdir}/pcnv/wigs") or die "$!";

	}
	
	if( @jobids > $cfg_hash{"NJOBS_LIMIT"} ){
	
		$jobIDString = join(',', @jobids);
	
		print "\n\tJob IDs: $jobIDString";
		
		$stdout = `echo "" | $cfg_hash{"QSUB"} -V -N MONITOR.TRAININGBAM2WIG -hold_jid $jobIDString -q ${queue}`;
		@tp = split(' ', $stdout);
		$holdJID = $tp[2];
		
		@jobids='';

	}
	 
}
close ( IN );

push(@jobids, $holdJID);
$jobIDString = join(',', @jobids);
$stdout = `echo "cd ${outdir}/pcnv; touch sortTrainingBam.complete" | $cfg_hash{"QSUB"} -V -N CHECKPOINT.TRAININGBAM2WIG -m ae -M $email -q ${queue} -hold_jid ${jobIDString}`;
@tp = split(' ', $stdout);
$holdJID = $tp[2];




# --------------------------------------------------------------------------------
# define stage-specific operations
# --------------------------------------------------------------------------------
sub stage1_pcnv{

	&printToLogNStdout( " (PatternCNV)");
	
	my $outdir_pcnv = "${outdir}/pcnv";
	$errors=0;

	# Remove breadcrumb files if present
	 unlink "${outdir_pcnv}/sortbam.complete";
	 unlink "${outdir_pcnv}/pcnv.complete";
	 unlink "${outdir_pcnv}/pcnv.failed";

	# create output dir # system call
	  `mkdir -p "${outdir_pcnv}"`;

	# create sample info file
	 open (OUT, ">${outdir}/configs/sample_info.txt") or die "\n\tERROR: Unable to create sample_info.txt file in ${outdir}/configs";
		# add header
 		  print OUT "sample.name\tsubject.ID\tsample.type\tbatch.ID\tBAM.file";
		  print OUT "\tsex" if $genderKnown;
		# add training set sample info
		  open(IN, "< $cfg_hash{TRAININGSET_BAMS}") or die $!;
		  while(<IN>){
			chomp;

			my $bamPath = $_;
			my $bamBaseName = basename($bamPath);
			my @tp = split('\.',$bamBaseName);
			my $sampleID = $tp[0];

			print OUT "\n${sampleID}\t${sampleID}\tGermline\t1\t${bamPath}";

			if($genderKnown){
				if( defined $sample_gender{$sampleID} ){
					print OUT "\t$sample_gender{$sampleID}";
				}
				else{
					print "\n\tError: Gender info for training sample $sampleID not specified.";
					$errors = $errors + 1;
				}
			}
			
		  }
		  close(IN);
		# add target sample info
		  foreach my $sampleID (@SAMPLEIDS){
			#print OUT "\n${sampleID}\t${sampleID}\tSomatic\t1\t${indir}/${sampleID}\.sorted\.bam";
			print OUT "\n${sampleID}\t${sampleID}\tSomatic\t1\t${indir}/${sampleID}\.bam";

			if($genderKnown){
				if( defined $sample_gender{$sampleID} ){
					print OUT "\t$sample_gender{$sampleID}";
				}
				else{
					print "\n\tError: Gender info for target sample $sampleID not specified.";
					$errors = $errors + 1;
				}
			}

		  }
	 close (OUT);

	# Check for errors
	 if($errors){
		die "\nThere were error(s) running stage 1 (PCNV).";
		return 2;
	 }
	 
	# write pcnv config file
	  open ( OUT, "> $outdir/configs/pcnv.config.txt" ) or die "\nSTAGE 1 ERROR: $!\n\n";
		print OUT "## PatternCNV";
		print OUT "\nOUTPUT_DIR=$outdir_pcnv";
		print OUT "\nSAMPLE_INFO=${outdir}/configs/sample_info.txt";
		print OUT "\nCAPTUREKIT_BED=$cfg_hash{CAPTUREKIT_BED}";
		print OUT "\nEXON_BED=$cfg_hash{EXON_BED}";
		print OUT "\nGENOME_SIZE=$cfg_hash{GENOME_SIZE}";
		print OUT "\nGENOME_REF=$cfg_hash{GENOME_REF}";
		print OUT "\n\n## Tools";
		print OUT "\nPATTERNCNV=$cfg_hash{PATTERNCNV}";
		print OUT "\nSAMTOOLS=$cfg_hash{SAMTOOLS}";
		print OUT "\nBEDTOOLS=$cfg_hash{BEDTOOLS}";
		print OUT "\nR=$cfg_hash{R}";
		print OUT "\n\n## SGE Parameters";		## may need modification for use on CGSL servers
		print OUT "\nEMAIL=${email}";
		print OUT "\nQUEUE=${queue}";				
		print OUT "\nQSUB_EXONKEY_MEMORY=$cfg_hash{QSUB_EXONKEY_MEMORY}";
		print OUT "\nQSUB_BAM2WIG_MEMORY=$cfg_hash{QSUB_BAM2WIG_MEMORY}";
		print OUT "\nQSUB_CALLCNVS_MEMORY=$cfg_hash{QSUB_CALLCNVS_MEMORY}";
	  close( OUT );
	 
	# launch pcnv and collect job ids
	  &printToLogNStdout("\n\tQsub PCNV jobs...");
	  my @stdout_pcnv = ();
	  if ($debug)	{
		@stdout_pcnv = `"$cfg_hash{PATTERNCNV}"/patternCNV_wrapper.sh -c $outdir/configs/pcnv.config.txt -b $cfg_hash{COVG_RES} -z $cfg_hash{BIN_SIZE} -m $cfg_hash{MIN_MAPPING_QUAL} -x $cfg_hash{EXTENSION_BUFFER} -d`;
	  } else {
		print "$cfg_hash{PATTERNCNV}/patternCNV_wrapper.sh -c $outdir/configs/pcnv.config.txt -b $cfg_hash{COVG_RES} -z $cfg_hash{BIN_SIZE} -m $cfg_hash{MIN_MAPPING_QUAL} -x $cfg_hash{EXTENSION_BUFFER} -d";
		@stdout_pcnv = `"$cfg_hash{PATTERNCNV}"/patternCNV_wrapper.sh -c $outdir/configs/pcnv.config.txt -b $cfg_hash{COVG_RES} -z $cfg_hash{BIN_SIZE} -m $cfg_hash{MIN_MAPPING_QUAL} -x $cfg_hash{EXTENSION_BUFFER}`;
	  }
	  &printToLogNStdout( "\n\tPCNV jobs launched: ");
	  my $startTime = time();
	  my $elapsedTime = 0;

	# extract job ids from pcnv stdout
	  @jobids = ();
	  foreach (@stdout_pcnv){
		if(/^Your job .+/){
			my @tp = split(' ', $_);
			push(@jobids, $tp[2]);
		}
	  }
	  push(@jobids, $holdJID);
	  $jobIDString = join(',', @jobids); 
	  &printToLogNStdout( $jobIDString);
	  
	  
	# create breadcrumb file upon job completion
		 $stdout = `echo "cd ${outdir}; touch PCNV.complete" | $cfg_hash{"QSUB"} -V -m ae -M $email -N CHECKPOINT.PATTERNCNV -q ${queue} -wd ${outdir}/logs/jobs -hold_jid $jobIDString`;
		 @tp = split(' ', $stdout);
		 $holdJID = $tp[2];
	  
	 
}




sub stage2_QC{

	&printToLogNStdout( " (QC)");

	my $outdir_qc = "${outdir}/qc"; `mkdir -p $outdir_qc`;

	# generate QC metrics - MAD diff cnvLg2Ratios, STDEV cnvLg2Ratios, Corr., Figure WRT training set

	 &printToLogNStdout( "\n\t* Generating STDEV, MAD lg2CNV matrices...");
	 $exec="$cfg_hash{R}/Rscript $cfg_hash{PIPELINE}/scripts/gen.QC.metrics.wrtTrainingSet.R ${outdir_qc} ${outdir}/pcnv/cnv-txt/Somatic_CNV_matrix.txt ${outdir}/pcnv/cnv-txt/Germline_CNV_matrix.txt ${outdir}/pcnv/cnv-txt/Somatic_RPKM_matrix.txt ${outdir}/pcnv/cnv-txt/Germline_RPKM_matrix.txt NA ${outdir}/pcnv/wigs $cfg_hash{MAX_MAD_LG2RATIO_DIFF} $cfg_hash{MAX_STDEV_LG2RATIO} $cfg_hash{MIN_CORR} $cfg_hash{MIN_COVG} $cfg_hash{MAX_NUM_DELS} $cfg_hash{MAX_NUM_DUPS} $cfg_hash{CUTOFF_LG2R_DUPS} $cfg_hash{CUTOFF_LG2R_DELS} $cfg_hash{CUTOFF_ZSCORE_DUPS} $cfg_hash{CUTOFF_ZSCORE_DELS} $indir $cfg_hash{CANDIDATE_TRAININGSET_BAMS}";
	 &printToLogNStdout("\n\t$exec");
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N GEN.QCTABLE -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid $holdJID`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];

}




sub stage3_tertiaryAnalyses{

	&printToLogNStdout( " (Tertiary Analyses)");

	my $outdir_tert = "${outdir}/tert"; 
	`mkdir -p $outdir_tert`;
	`mkdir -p ${outdir}/VCF`;

	
	# generate pVals for each segment
	 &printToLogNStdout( "\n\n\t* Generating p-vals for each segment...");
	 $exec="$cfg_hash{R}/Rscript $cfg_hash{PIPELINE}/scripts/gen.seg.pVal.matrix.R ${outdir}/pcnv/cnv-txt $cfg_hash{GAINS_HIGH_TRANCHE_MIN_PVAL} $cfg_hash{GAINS_HIGH_TRANCHE_MAX_PVAL} $cfg_hash{GAINS_MEDIUM_TRANCHE_MIN_PVAL} $cfg_hash{GAINS_MEDIUM_TRANCHE_MAX_PVAL} $cfg_hash{GAINS_LOW_TRANCHE_MIN_PVAL} $cfg_hash{GAINS_LOW_TRANCHE_MAX_PVAL} $cfg_hash{LOSSES_HIGH_TRANCHE_MIN_PVAL} $cfg_hash{LOSSES_HIGH_TRANCHE_MAX_PVAL} $cfg_hash{LOSSES_MEDIUM_TRANCHE_MIN_PVAL} $cfg_hash{LOSSES_MEDIUM_TRANCHE_MAX_PVAL} $cfg_hash{LOSSES_LOW_TRANCHE_MIN_PVAL} $cfg_hash{LOSSES_LOW_TRANCHE_MAX_PVAL}";
	 &printToLogNStdout( "\n\t$exec");
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N GEN.PVALS -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid ${holdJID}`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];
	
	
	# generate exon-level VCF
	 &printToLogNStdout( "\n\t* Generating VCF file (for exons)...");
	 $exec="$cfg_hash{R}/Rscript $cfg_hash{PIPELINE}/scripts/gen.VCF_exon.R $project $email $cfg_hash{GENOME_REF_NAME} ${outdir}/pcnv/cnv-txt/Somatic_CNV_matrix.txt ${outdir}/pcnv/configs/exon_key.txt ${outdir}/qc/qc.metrics.csv ${sampleMetaData} $cfg_hash{CUTOFF_LG2R_DELS} $cfg_hash{CUTOFF_LG2R_DUPS} ${outdir}/VCF";
	 &printToLogNStdout( "\n\t$exec");
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N GEN.EXONVCF -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid ${holdJID}`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];

	 
	# generate segment-level VCF
	 &printToLogNStdout( "\n\t* Generating VCF file (for CBS segments)...\n");
	 $exec="$cfg_hash{R}/Rscript $cfg_hash{PIPELINE}/scripts/gen.VCF_segment.R $project $email $cfg_hash{GENOME_REF_NAME} ${outdir}/pcnv/cnv-txt/Somatic_CNV_matrix.txt ${outdir}/pcnv/configs/exon_key.txt ${outdir}/qc/qc.metrics.csv ${sampleMetaData} $cfg_hash{CUTOFF_LG2R_DELS} $cfg_hash{CUTOFF_LG2R_DUPS} ${outdir}/VCF ${outdir}/pcnv/cnv-txt";
	 &printToLogNStdout( "\n\t$exec");
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N GEN.SEGMENTVCF -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid ${holdJID}`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];
	 
	 
	# merge segment level VCFs and reformat merged file for upload to VCF-Miner 
	 &printToLogNStdout( "\n\t* Merging/reformatting segment VCF files (for CBS segments)...\n");
	 $exec="perl /data2/bsi/staff_analysis/m139467/vcf_scripts/mergeSegmentVCFs.pl \$(ls -m ${outdir}/VCF/*.segment.vcf | tr -d ' ' | tr -d '\\n' | tr -d '\\r') ${outdir}/VCF/${project}_merged.segment.vcf;";
	 $exec.="perl /data2/bsi/staff_analysis/m141962/chenWang/SL2/vcfFormat2Info/scripts/vcf_format2info.pl ${outdir}/VCF/${project}_merged.segment.vcf ${outdir}/VCF/${project}_merged.segment.format2info.vcf";
	 &printToLogNStdout( "\n\t$exec");
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N GEN.MERGED.SEGMENT.FORMAT2INFOVCF -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid ${holdJID}`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];
	 
	 
	# # reformat segment-level VCFs for loading into VCF-Miner
	 # &printToLogNStdout( "\n\t* Reformatting segment VCFs for VCF-Miner...\n");
	 # $exec="find ${outdir}/VCF -name '*.segment.vcf' | while read VCF; do outputFile=\"\${VCF/segment.vcf/segment.format2info.vcf}\"; /data2/bsi/staff_analysis/m141962/chenWang/SL2/SL2_pipeline/scripts/vcf_format2info.pl \$VCF \$outputFile; done";
	 # &printToLogNStdout( "\n\t$exec");
	 # $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N REFORMAT.SEGMENTVCF -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid ${holdJID}`;
	 # &printToLogNStdout("\n\t$stdout");
	 # @tp = split(' ', $stdout);
	 # $holdJID = $tp[2];


}




sub stage4_prepareDeliveryFolder{
	
	 &printToLogNStdout( " (Prepare Delivery Folder)\n" );
	 
	 my $outdir_delivery = "${outdir}/cnvDelivery"; 
	 `mkdir -p $outdir_delivery`;
	 
	 $exec="$cfg_hash{R}/Rscript $cfg_hash{PIPELINE}/scripts/create_cnv_txt.R ${outdir}/pcnv/cnv-txt $outdir_delivery ;";
	 $exec.="cp -r ${outdir}/VCF ${outdir_delivery};";
	 $exec.="cp ${outdir}/pcnv/configs/exon_key.txt ${outdir_delivery};";
	 $exec.="cp ${outdir}/pcnv/cnv-txt/*_CNV_seg_info.txt ${outdir_delivery};";
	 $exec.="cp ${indir}/../docs/config/analysis.ped ${outdir_delivery};";
	 $exec.="cp -r ${outdir}/pcnv/wigs ${outdir_delivery};";

	 &printToLogNStdout("\n\t$exec");
	 
	 $stdout=`echo "$exec" | $cfg_hash{QSUB} -m abe -M $email -N PREPARE.DELIVERYFOLDER -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid $holdJID`;
	 &printToLogNStdout("\n\t$stdout");
	 @tp = split(' ', $stdout);
	 $holdJID = $tp[2];
 
 
}




my %stage_operations = ( 1 => \&stage1_pcnv,
			 2 => \&stage2_QC,
			 3 => \&stage3_tertiaryAnalyses,
			 4 => \&stage4_prepareDeliveryFolder );

			  
			  
			  
# --------------------------------------------------------------------------------
# Execute user-specified stages
# --------------------------------------------------------------------------------
foreach my $stageID ( sort { $a <=> $b } @stages ){
	$time = localtime;
	&printToLogNStdout( "\n\n\n[$time] BEGIN STAGE $stageID");
	&{${stage_operations{$stageID}}};
	
}

`echo "cd ${outdir}; touch pipeline.complete" | $cfg_hash{QSUB} -m ae -M $email -N CHECKPOINT.PIPELINE -V -q ${queue} -wd ${outdir}/logs/jobs -hold_jid $holdJID`;
&printToLogNStdout( "\n\n\n[$time] ALL JOBS SUBMITTED ---");
print "\n\n\n";




# --------------------------------------------------------------------------------
# HELPER FUNCTIONS
# --------------------------------------------------------------------------------

sub usage { print "\n  Usage: $SCRIPTNAME";
            print "\n\t\t -id indir             # REQUIRED path/to/samplefolder";
			print "\n\t\t -od outdir            # REQUIRED path/to/outputfolder";
            print "\n\t\t -pj project           # REQUIRED project name";
            print "\n\t\t -mc mainconfig        # REQUIRED config file for PatternCNV";
			print "\n\t\t -pc captureKitConfig  # REQUIRED capture-kit specific config file";
            print "\n\t\t -sd sampleMetaData    # OPTIONAL sample metadata file";
			print "\n\t\t  -q queue             # queue; DEFAULT sandbox.q";
			print "\n\t\t  -e email             # email; DEFAULT <you>\@mayo.edu";
			print "\n\t\t  -s stage             # (a) all, (1) pcnv, (2) qc, (3) tertiaryAnalyses, (4) prepareDeliveryFolder; DEFAULT: a";
			print "\n\t\t -dg debug mode        # (0) off, (1) on; DEFAULT 0";
			print "\n\t\t -rw regenerate wigs   # (0) off, (1) on; DEFAULT 0";
			print "\n\t\t  -h display usage\n\n"   ;
}

sub quit{
	my $exitCode = shift;
	&printToLogNStdout( "\n\n  Pipeline terminated with errors! (exit code: $exitCode)\n\n");
	exit $exitCode;
}

sub progressAnimation{
	#todo
}

sub lastName {
	(my $sub_input) = @_; 
	my @sub_inline = split(/\//, $sub_input); 
	my $sub_fname = $sub_inline[scalar(@sub_inline)-1]; 
	return $sub_fname;
}

sub sleepforjob {
  my ($sub_max, $sub_jobname, $sub_sleep ) = @_;
  my $sub_w = `whoami`; chomp ( $sub_w );
  my $sub_jobline = `$cfg_hash{QSTAT} | grep $sub_w | grep $sub_jobname`;

  my @sub_line = split ( /\n/, $sub_jobline );
  my $sub_job = scalar(@sub_line);
  if ( $sub_job > $sub_max ) { sleep $sub_sleep; }
  return $sub_job;
}

sub printToLogNStdout {
	my $text = shift;
	print $text;
	print LOG $text;

}


__END__



