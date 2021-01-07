#!/usr/bin/env perl

# Extract Coverage from bam files
use strict;
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long;
use Config;

my $dirname = dirname(__FILE__);

my $samtools;
my $bedtools;
my $mosdepth;
my $bedGraphToBigWig;
my $genomeVersion;

my $Rscript = `which Rscript`;
chomp $Rscript;
if (!$Rscript){
	print " ERROR: Rscript was not found on PATH\n";
} 

# Setting binaries for either mac (darwin) or Linux systems
if ($Config{osname} =~/darwin/) {

	$bedGraphToBigWig =( -x "$dirname/bin/darwin/bedGraphToBigWig" )
	? "$dirname/bin/darwin/bedGraphToBigWig"
	:  die " ERROR: Unable to execute bedGraphToBigWig"; 

	$mosdepth = ( -x "$dirname/bin/darwin/mosdepth" )
	? "$dirname/bin/darwin/mosdepth"
	:  die " ERROR: Unable to execute mosdepth";

	$bedtools = ( -x "$dirname/bin/darwin/bedtools" )
	? "$dirname/bin/darwin/bedtools"
	:  die " ERROR: Unable to execute bedtools";

	$samtools = ( -x "$dirname/bin/darwin/samtools" )
	? "$dirname/bin/darwin/samtools"
	:  die " ERROR: Unable to execute samtools";
}
else {

	$bedGraphToBigWig =( -x "$dirname/bin/linux/bedGraphToBigWig" )
	? "$dirname/bin/linux/bedGraphToBigWig"
	:  die " ERROR: Unable to execute bedGraphToBigWig"; 

	$mosdepth = ( -x "$dirname/bin/linux/mosdepth" )
	? "$dirname/bin/linux/mosdepth"
	:  die " ERROR: Unable to execute mosdepth";

	$bedtools = ( -x "$dirname/bin/linux/bedtools" )
	? "$dirname/bin/linux/bedtools"
	:  die " ERROR: Unable to execute bedtools";

	$samtools = ( -x "$dirname/bin/linux/samtools" )
	? "$dirname/bin/linux/samtools"
	:  die " ERROR: Unable to execute samtools";
}

 my $bamDir;
 my $bedDir;
 my $minCov = 3;
 my $threads = 4;
 my $chromFile;
 my $normFactor = 10e7;
 my $doHistogram = "";

 Help () if (@ARGV < 1 or !GetOptions(
	'input|i=s'=>\$bamDir,
	'outdir|o=s'=>\$bedDir,
	'mincov|m=i'=>\$minCov,
	'normfactor|n=i'=>\$normFactor, 
	'threads|t=i'=>\$threads,
	'reference|r=s'=>\$genomeVersion
	)
 );

 if (!$genomeVersion) {
	 print " ERROR: Missing genome version (hg19 or hg38)\n";
	 Help();
 }
 if ($genomeVersion ne 'hg19' && $genomeVersion ne 'hg38') {
	 print " ERROR: Missing genome version (hg19 or hg38)\n";
	 Help();
 }
 if ($genomeVersion eq 'hg19') {
	$chromFile = "$dirname/genome/hg19.chrom.sizes";
 }
 if ($genomeVersion eq 'hg38') {
	$chromFile = "$dirname/genome/hg38.chrom.sizes";
 }

 if (!-e $bamDir) {
	print " ERROR: $bamDir does not exist\n"; exit;
 }
 mkdir $bedDir;

 my @bams = glob ("$bamDir/*bam");
 if (!@bams) {
     print " ERROR: no BAMs found on $bamDir directory\n"; exit;
 }

 my $pm = Parallel::ForkManager->new($threads);

 my $localtime = localtime();        
 print "$localtime - INFO: Using a total number of $threads CPUs\n";
 print "$localtime - INFO: Extracting regions with at least $minCov reads\n";

 # Creating BAM indices if not available
 foreach my $bam ( @bams) {
    my $bai = $bam . ".bai";
    if (!-e $bai) {
	    my $localtime = localtime();        
		print "$localtime - INFO: Indexing $bam\n";

        my $cmd = "$samtools index $bam";
		system $cmd;
    }
 }

 # Coverage extraction
 foreach my $bam ( @bams ) {

    my $localtime = localtime();
    print "$localtime - INFO: Processing $bam\n";

    my $name = basename ($bam);
    $name =~s/.bam//;
	my $perBaseCov = "$bedDir/$name.per-base.bed.gz";

    # Extracting per base coverage with mosdepth
	my $cmd = "$mosdepth --fast-mode --mapq 50 -t $threads $bedDir/$name $bam";
	system $cmd if !-e $perBaseCov;

	my $plainPerBaseCov = $perBaseCov;
	$plainPerBaseCov =~s/.bed.gz/.bed/;

    $localtime = localtime();
    print "$localtime - INFO: Sorting per base depth file for sample $name\n";
	$cmd = "gunzip -c $perBaseCov | LC_COLLATE=C sort -k1,1 -k2,2n > $plainPerBaseCov";
	system $cmd if !-e $plainPerBaseCov;

	if (-e $perBaseCov && !-z $perBaseCov) {

	    # Fast count of total reads. 
	    # Note that since we have filtered out secondary/suplementary reads this fast method is similar than doing "view -c"	
		my $counts = `$samtools idxstats $bam | awk -F '\t\' \'{s+=\$3+\$4}END{print s}'`;
		chomp $counts;

		# Getting regions with enough coverage
		$cmd = "zcat $perBaseCov | awk '{ end=\$3; if (\$4 >= $minCov) {  print \$1\"\t\"\$2\"\t\"end } }'".
		" | $bedtools merge -i stdin > $bedDir/$name.COVERAGE.bed";
		system $cmd if !-e "$bedDir/$name.COVERAGE.bed";

		# Creating bigwig file using UCSC utilities
		$localtime = localtime();
		print "$localtime - INFO: Creating bigWig for sample $name\n";
		bedGraph2BigWig($plainPerBaseCov, $normFactor, $counts);
		unlink($plainPerBaseCov);
	}
	else {
		print " ERROR: $bedDir/$name.per-base.bed.gz was not created\n";
		exit;
	}
 }
##################################
sub doHistogram {
	my $covFiles = shift;
	my $regions  = shift;
} 

##################################
sub bedGraph2BigWig {

	my $bedGraph   = shift;
	my $normFactor = shift;
	my $counts     = shift;

	my $scale = $normFactor/$counts;
	my $bigWig   = $bedGraph;
	$bigWig =~s/.per-base//;
	$bigWig =~s/.bed/.bw/;

	my $tmpBedGraph = $bedGraph;
	$tmpBedGraph =~s/.bed/.normalized.bed/;
	open (IN, "<", $bedGraph) or die " ERROR: Unable to open $bedGraph\n";
	open (OUT, ">", $tmpBedGraph) or die " ERROR: Unable to open $tmpBedGraph\n";
	while (my $line=<IN>){
		chomp $line;
		my @tmp = split ("\t", $line);
		my $counts = $tmp[-1];
		if ($counts > 0){
			$tmp[-1] = $tmp[-1]*$scale;  
		}
		my $newline = join("\t",@tmp);
		print OUT "$line\n";
	} 
	close IN;
	close OUT;

	my $cmd = "$bedGraphToBigWig $tmpBedGraph $chromFile $bigWig";
	system $cmd if !-e $bigWig;

	$cmd = "gzip -f $tmpBedGraph";
	system $cmd if !-e $tmpBedGraph . ".gz";
}
##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -i,--input       STRING   Input directory with BAM files
 -o,--outdir      STRING   Output directory
 -m,--mincov      INT      Minimum coverage required at each sequenced position (default=3)
 -n,--normfactor  INT      Normalization factor (default=10e7)
 -r,--reference   STRING   Genome versions (choose: hg19, hg38, default=hg19)
 -t,--threads     INT      Number of CPU cores (default=4)\n\n";
 exit;
}

