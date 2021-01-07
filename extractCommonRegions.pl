#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);
use List::Util qw(first max);
use Getopt::Long;
use Config;

my $dirname = dirname(__FILE__);
my $inputDir;
my $THREADS = 4;
my $devNull = ">/dev/null 2>&1";

 Help () if (@ARGV < 2 or !GetOptions(
	'input|i=s'=>\$inputDir,
	'threads|t=i'=>\$THREADS
	)
 );

if (!$inputDir) {
	print " ERROR: No input dir was specified\n";
	exit;
}
my @bedCovs = glob ("$inputDir/*COVERAGE.bed");
if (!@bedCovs) {
	print " ERROR: No input *COVERAGE.bed files were found on $inputDir directory\n";
	exit;
}
my $pm = Parallel::ForkManager->new($THREADS);
my $bedtools;

# Setting binaries for either mac (darwin) or Linux systems
if ($Config{osname} =~/darwin/) {
	$bedtools = ( -x "$dirname/bin/darwin/bedtools" )
	? "$dirname/bin/darwin/bedtools"
	:  die " ERROR: Unable to execute bedtools";
}
else {
	$bedtools = ( -x "$dirname/bin/linux/bedtools" )
	? "$dirname/bin/linux/bedtools"
	:  die " ERROR: Unable to execute bedtools";
}

our %overlaps = ();
my  @arrayOlaps;
my $localTime = localtime();

my $tmpOlaps = "$inputDir/tmpOlaps.txt";
open (TMP, ">", $tmpOlaps) || die " ERROR: Unable to open $tmpOlaps\n";

foreach my $bed1 (@bedCovs) {

	my $name1 = basename ($bed1);
	$name1=~s/.rmdup.bam.COVERAGE.bed//;
	$name1=~s/.COVERAGE.bed//;

	foreach my $bed2 (@bedCovs) {
	   $pm->start and next;

		my $name2 = basename ($bed2);
		$name2=~s/.rmdup.bam.COVERAGE.bed//;
	    $name2=~s/.COVERAGE.bed//;

		if (!exists $overlaps{$name1}{$name2} && !exists $overlaps{$name2}{$name1} && $name1 ne $name2) { 

            $localTime = localtime();
            print "$localTime - INFO: Intersecting $name1 and $name2\n";

			my $olap  = `$bedtools intersect -a $bed1 -b $bed2 -wao -sorted -nonamecheck | awk '{x += \$7} END {print x}'`;
			chomp $olap;
			push @arrayOlaps, $olap;

			my $cmd = "$bedtools intersect -a $bed1 -b $bed2  -sorted -nonamecheck > $inputDir/$name1.$name2.shared.bed";
			system $cmd;

			$overlaps{$name1}{$name2} = $olap;
			$overlaps{$name2}{$name1} = $olap;

			print TMP "$name1\t$name2\t$overlaps{$name1}{$name2}\n";
		}	
		$pm->finish;
	}
	$pm->wait_all_children;
}

open (TMP, "<", $tmpOlaps) || die " ERROR: Unable to open $tmpOlaps\n";
while (my $line =<TMP>) {
	chomp $line;
	my @tmp = split (/\t/, $line);
	my $name1 = $tmp[0];
	my $name2 = $tmp[1];
	my $olap = $tmp[2];

	$overlaps{$name1}{$name2} = $olap;
    $overlaps{$name2}{$name1} = $olap;

    push @arrayOlaps, $olap;
}
close TMP;
unlink($tmpOlaps);

##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -i,--input     STRING   Input directory with files ending with *COVERAGE.bed
 -t,--threads   INT      Number of CPU cores (default = 4)\n\n";
	exit;
}