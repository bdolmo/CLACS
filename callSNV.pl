#!/usr/bin/env perl

use strict;
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long;
use Config;

my $dirname = dirname(__FILE__);

my $bwa;
my $samtools;
my $freebayes;
my $gatk;
my $devNull = ">/dev/null 2>&1";

# Setting binaries for either mac (darwin) or Linux systems
if ($Config{osname} =~/darwin/) {
	$samtools = ( -x "$dirname/bin/darwin/samtools" )
	? "$dirname/bin/darwin/samtools"
	:  die " ERROR: Unable to execute samtools";

    $freebayes = ( -x "$dirname/bin/darwin/freebayes" )
	? "$dirname/bin/darwin/freebayes"
	:  die " ERROR: Unable to execute freebayes";
}
else {
	$samtools = ( -x "$dirname/bin/linux/samtools" )
	? "$dirname/bin/linux/samtools"
	:  die " ERROR: Unable to execute samtools";

    $freebayes = ( -x "$dirname/bin/linux/freebayes" )
	? "$dirname/bin/linux/freebayes"
	:  die " ERROR: Unable to execute freebayes";
}

# $gatk = "$dirname/bin/linux/GenomeAnalysisTK.jar";
# if (!$gatk) {
#     print " ERROR: Unable to execute GATK\n";
#    exit;
# }

my $bamDir;
my $vcfDir;
my $reference;
my $threads = 4;

Help () if (@ARGV < 1 or !GetOptions(
	'input|i=s'=>\$bamDir,
	'outdir|i=s'=>\$vcfDir,
	'genome|g=s'=>\$reference,
	'threads|t=i'=>\$threads
	)
);

my $localTime = localtime();

if (!-e $bamDir) {
    print " ERROR: input bam dir (-i) $bamDir does not exist\n";
    exit;
}
if (!-e $reference) {
    print " ERROR: reference genome (-r) $reference does not exist\n";
    exit; 
}
my $refDict = $reference;
$refDict =~s/.fasta/.dict/;

if (!-e $refDict){
    print "$localTime - INFO: Creating dict file from reference $reference\n";
    my $cmd = "$samtools dict $reference -o $refDict";
    system $cmd;
}

mkdir $vcfDir;

my $pm = Parallel::ForkManager->new($threads);
my @bams = glob ("$bamDir/*.bam");

if (!@bams) {
    print " ERROR: It seems that there are no bam files (prefix = .bam) on $bamDir directory\n";
    exit;
}
foreach my $bam (@bams) {
    $pm->start and next;
    $localTime = localtime();

    my $name = basename($bam);
    $name =~s/.nodup.bam//;
    $name =~s/.bam//;

    # print "$localTime - INFO: Variant call (GATK HC) over sample $name\n";
    print "$localTime - INFO: Calling variants (FreeBayes) on sample $name\n";

    # my $cmd = "java -jar $gatk -T HaplotypeCaller -I $bam -o $vcfDir/$name.vcf -R $reference $devNull";
    # system $cmd if !-e "$vcfDir/$name.vcf";
    #print "$cmd\n";
    my $cmd = "$freebayes -f $reference $bam > $vcfDir/$name.raw.vcf";
    system $cmd if !-e "$vcfDir/$name.raw.vcf";
    $pm->finish;
}
$pm->wait_all_children;

# Now filter vcfs
my @vcfs = grep ($_!~/filtered/, glob ("$vcfDir/*.vcf"));
foreach my $vcf (@vcfs){
    # filterGatkVCF($vcf);
    my $filteredVCF = filterFreeBayesVCF("$vcf");
    unlink("$vcf");
}  

##################################

sub filterGatkVCF {

    my $rawVCF = shift;
    my $filteredVCF = $rawVCF;
    $filteredVCF =~s/.vcf$/.filtered.vcf/;   
    open (IN, "<", $rawVCF);
    open (OUT, ">", $filteredVCF);
    while (my $line=<IN>){
        chomp $line;
        my@tmp = split ("\t", $line);
        if ($line =~/^#/){
            print OUT "$line\n";
            next;
        }  
        else{
            my $ref = $tmp[3]; 
            my $alt = $tmp[4];
            my $qual= $tmp[5];
            my $filter = $tmp[6]; 

            # Skipping indels 
            if (length($ref) > 1 or length($alt) >1){
                next;
            } 

            # filtering by DP    
            my@info = split (";", $tmp[7] );
            my ($DP) = grep ($_=~/^DP=/,@info );
            $DP=~s/DP=// if $DP;

            # filtering by Map quality
            my ($MQ) = grep ($_=~/^MQ=/,@info );
            $MQ=~s/MQ=// if $MQ;

            # filtering by strand bias odds ratio
            my ($SOR) = grep ($_=~/^SOR=/,@info );
            $SOR=~s/SOR=// if $SOR;

            # filtering by strand bias odds ratio
            my ($QD) = grep ($_=~/^QD=/,@info );
            $QD=~s/QD=// if $QD;

            my $MQRankSum = grep ($_=~/^MQRankSum=/, @info );
            $MQRankSum=~s/MQRankSum=// if $MQRankSum;

            my $ReadPosRankSum = grep ($_=~/^ReadPosRankSum=/, @info );
            $ReadPosRankSum=~s/ReadPosRankSum=// if $ReadPosRankSum;

            my @format_tag = split(":", $tmp[8]);
            my @format  = split(":", $tmp[9]);
            my $AD = $format[1];
            my @tmp_reads = split(",", $AD);
            my $alt_reads = $tmp_reads[1]; 

            if ($MQRankSum){
                if ($MQRankSum < -12.5){
                    $filter = "lowQual";
                } 
            } 
            if ($ReadPosRankSum){
                if ($ReadPosRankSum < -8){
                    $filter = "lowQual";
                } 
            } 
            if ($DP < 5 or $MQ < 40 or $SOR > 3 or $QD < 2 or $alt_reads < 2){
                $filter = "lowQual";
            }  
            else{
                $filter = "PASS";
            }
            $tmp[6] = $filter;
            $line = join("\t", @tmp); 
        } 
        print OUT "$line\n"; 
    } 
    close IN;
    close OUT;

    unlink($rawVCF);
    rename $filteredVCF, $rawVCF;
} 

##################################
sub filterFreeBayesVCF {

    my $rawVCF = shift;
    my $filteredVCF = $rawVCF;
    $filteredVCF =~s/.raw.vcf/.vcf/;   
    open (IN, "<", $rawVCF);
    open (OUT, ">", $filteredVCF);
    while (my $line=<IN>){
        chomp $line;
        my@tmp = split ("\t", $line);
        if ($line =~/^#/){
            print OUT "$line\n";
            next;
        }  
        else{
            my $ref = $tmp[3]; 
            my $alt = $tmp[4];
            my $qual= $tmp[5];
            next if $qual < 20;

            # Skipping indels 
            if (length($ref) > 1 or length($alt) >1){
                next;
            } 

            # Now filtering by DP    
            my@info = split (";", $tmp[7] );
            my ($DP) = grep ($_=~/^DP=/,@info );
            $DP=~s/DP=//; 
            #next if $DP < 5;
        } 
        print OUT "$line\n"; 
    } 

    close IN;
    close OUT;
} 

##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -i,--input    STRING   Input directory with bam files
 -o,--outdir   STRING   Output directory
 -g,--genome   STRING   Reference genome in FASTA format
 -t,--threads  INT      Number of CPU cores (default = 4)\n\n";
	exit;
}