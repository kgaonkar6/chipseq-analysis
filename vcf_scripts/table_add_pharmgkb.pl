
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


my($input, $sample, $filter, $help,$outvcf,$inrsID);
GetOptions("in|i=s"     => \$inrsID,
           "filter|f=s"    => \$filter,
           "outvcf|o=s"    => \$outvcf)


`for ID in $(cat $inrsID ); do grep -w $ID variants.pharmgkb.tsv|cut -f4|cut -d "," -f1|cut -d ":" -f2; done >${inrsID}.filtered`
open FILT_rsID ,${inrsID}.filtered, or die "opening ${inrsID}.filtered: $!\n"
open FILT_rsIDJOB, ${inrsID}.filtered.job,or die "opening ${inrsID}.filtered.job: $!\n"

while (my $ID =<FILT_rsID> ){
print FILT_rsIDJOB "echo". $ID ."|bior_lookup -d /data5/bsi/catalogs/user/v1/pharmacogenetics/pharmgkb/2013_07_03/variants.pharmgkb.tsv.bgz -p RSID |bior_drill -p DRUGS -p DISEASES";
}
close FILT_rsIDJOB;

`paste <(for id in $(source ${inrsID}.filtered.job |grep -v "#" |cut -f1);do awk -v ID=$id '{if( $21 == ID) print $0}' test.xls ;done ) <(source ${inrsID}.filtered.job |grep -v "#" |cut -f2-) > ${inrsID}.filtered.xls`

