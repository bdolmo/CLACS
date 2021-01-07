#!/usr/bin/env perl

use strict;
use File::Basename;
use Statistics::Basic qw(:all);
use Parallel::ForkManager;
use Getopt::Long;
use Config;

my $dirname = dirname(__FILE__);

 my $THREADS = 4;
 my $pm = Parallel::ForkManager->new($THREADS);

 my $vcfDir;
 my $bedDir;
 my $outDir;
 my $genomeVersion = "hg19";

 Help () if (@ARGV < 4 or !GetOptions(
	'vcf_dir|v=s'=>\$vcfDir,
	'bed_dir|b=s'=>\$bedDir,
	'reference|r=s'=>\$genomeVersion,
	'outdir|o=s'=>\$outDir,
	'threads|t=i'=>\$THREADS
	)
 );

 my @vcfs = glob ("$vcfDir/*.vcf");
 if (!@vcfs) {
     print " ERROR: no VCFs found on $vcfDir directory\n"; exit;
 }

 my @BEDS = glob ("$bedDir/*shared.bed");
 if (!@BEDS) {
     print " ERROR: no VCFs found on $bedDir directory\n"; exit;
 }

 mkdir $outDir;

 my %hash = ();
 my $samtools;
 my $bedtools;
 my %corr = ();
 my $blacklist;

 # Setting binaries for either mac (darwin) or Linux systems
 if ($Config{osname} =~/darwin/) {
	$bedtools = ( -x "$dirname/bin/darwin/bedtools" )
	? "$dirname/bin/darwin/bedtools"
	:  die " ERROR: Unable to execute bedtools";

	$samtools = ( -x "$dirname/bin/darwin/samtools" )
	? "$dirname/bin/darwin/samtools"
	:  die " ERROR: Unable to execute samtools";
 }
 else {
	$bedtools = ( -x "$dirname/bin/linux/bedtools" )
	? "$dirname/bin/linux/bedtools"
	:  die " ERROR: Unable to execute bedtools";

	$samtools = ( -x "$dirname/bin/linux/samtools" )
	? "$dirname/bin/linux/samtools"
	:  die " ERROR: Unable to execute samtools";
 }
 if ($genomeVersion eq 'hg19'){ 
    $blacklist = "$dirname/genome/ceph18.b37.exclude.2014-01-15.chr.bed";
 }
 else {
    $blacklist = "$dirname/genome/ceph18.b38.exclude.2014-01-15.chr.bed";
 }

if (!$genomeVersion) {
	print " ERROR: Missing genome version (hg19 or hg38)\n";
	Help();
}
if ($genomeVersion ne 'hg19' && $genomeVersion ne 'hg38') {
	print " ERROR: Incorrect genome version. Please, introduce hg19 or hg38\n";
	Help();
}

 my %samples = ();
 my %Genotypes = (

	"A/A" => 1,
	"A/T" => 2,
	"A/G" => 3,
	"A/C" => 4,

	"C/A" => 5,
	"C/C" => 6,
	"C/T" => 7,
	"C/G" => 8,

	"T/A" => 9,
	"T/C" => 10,
	"T/T" => 11,
	"T/G" => 12,

	"G/A" => 13,
	"G/C" => 14,
	"G/T" => 15,
	"G/G" => 16
);

my @samples = ();

foreach my $vcf1 ( @vcfs ) {
    my $name1 = basename ($vcf1);
    $name1 =~s/.vcf//;
    $name1 =~s/Aligned.sortedByCoord.out.COVERAGE.vcf//;
    push @samples, $name1;
}
my %seen = ();
open (TMP, ">", "$outDir/corr.tmp.txt") || die " ERROR: Unable to create $outDir/corr.tmp.txt\n";
foreach my $vcf1 ( @vcfs ) {

    my $name1 = basename ($vcf1);
    $name1 =~s/.vcf//;

    $samples{$name1}++;
    my $counter = 1;

    foreach my $vcf2 ( @vcfs ) {
        #$pm->start and next;
        $seen{$vcf1}{$vcf2}++;
        $seen{$vcf2}{$vcf1}++;

        next if $seen{$vcf1}{$vcf2} > 1;
        next if $seen{$vcf2}{$vcf1} > 1;

 	    $counter++;
        my $name2 = basename ($vcf2);
        $name2 =~s/.vcf//;
        if ($name1 eq $name2) {
            print TMP "$name1\t$name2\t.\t1\n";
        }
        else {
            my ($bed_candidate) = grep ($_ =~/$name1.$name2.shared.bed/, @BEDS);
            if (!$bed_candidate) {
                print " ERROR not found $name1 $name2 bed candidate\n";
            }
            else {

                my $localTime = localtime();
                print "$localTime - INFO: Comparing $name1-$name2 SNV profiles\n";

                my $cmd = "$bedtools intersect -a $vcf1 -b $blacklist -v -wa -header" .
                   " | $bedtools intersect -a stdin -b $bed_candidate -wa -header > $outDir/$name1.$counter.tmp.vcf";
                system $cmd if !-e "$outDir/$name1.$counter.tmp.vcf";

                $cmd = "$bedtools intersect -a $vcf2 -b $blacklist -v -wa -header" .
                   " | $bedtools intersect -a stdin -b $bed_candidate -wa -header > $outDir/$name2.$counter.tmp.vcf";
                system $cmd if !-e "$outDir/$name2.$counter.tmp.vcf";

                my $total_kb    = countKiloBases($bed_candidate);
                my $correlation = compareCoordinates($bed_candidate, 
                    "$outDir/$name1.$counter.tmp.vcf", "$outDir/$name2.$counter.tmp.vcf");
                print TMP "$name1\t$name2\t$total_kb\t$correlation\n";

                $localTime = localtime();
                print "$localTime - INFO: Correlation ($correlation) on a total of $total_kb\n";

                unlink("$outDir/$name1.$counter.tmp.vcf");
                unlink("$outDir/$name2.$counter.tmp.vcf");
            }
        } 
	    #$pm->finish();
    }
    #$pm->wait_all_children;
}
close TMP;

open (IN, "<", "$outDir/corr.tmp.txt") || die " ERROR: Unable to open $outDir/corr.tmp.txt\n";
while (my $line=<IN>) {
    chomp $line;
    my ($sample1, $sample2, $kb, $correlation) = split(/\t/, $line);
    $hash{$sample1}{$sample2} = $correlation;
    $hash{$sample2}{$sample1} = $correlation;
}
close IN;
unlink("$outDir/corr.tmp.txt");

open (OUT, ">", "$outDir/sample_correlations.tsv") 
    || die " ERROR: Unable to open $outDir/sample_correlations.tsv\n";

foreach my $sample1 ( sort keys %hash ) {
    print OUT "\t$sample1";
}
print OUT "\n";
foreach my $sample1 ( sort keys %hash ) {
    print OUT "$sample1";
    foreach my $sample2 ( sort keys %hash ) {
        print OUT "\t$hash{$sample1}{$sample2}";
    }
    print OUT "\n";
}
close OUT;
#unlink("$outDir/sample_correlations.tsv");

##############################
sub countKiloBases {
    my $bed = shift;
    open (IN, "<", $bed) || die " ERROR: Unable to open $bed\n";
    my $sum = 0;
    while (my $line=<IN>) {
        chomp $line;
        my @tmp = split (/\t/, $line);
        my $size = $tmp[2]-$tmp[1];
        $sum+=$size;
    }
    return $sum/1000 . " kb";
}

##############################
sub compareCoordinates {

    my $regions = shift;
    my $vcf_A   = shift;
    my $vcf_B   = shift;
    my %vcf = ();

    my @array_A = ();
    my @array_B = ();
    open (A, "<", $vcf_A) || die "ERROR: Unable to open $vcf_A\n";
    while (my $line=<A>) {
        chomp $line;
        next if $line =~/^#/;  
        my @tmp = split (/\t/, $line);
        my @info = split (";", $tmp[7]);
        next if $tmp[0] =~/_/;
        my $p = $tmp[1]+1;
	    #my $var = "$tmp[3]/$tmp[4]";
	    #$vcf{"$tmp[0]\t$tmp[1]\t$p"}{A} = $Genotypes{$var};
	    $vcf{"$tmp[0]\t$tmp[1]\t$p"}{A} = 1;		
    }
    close A;

    open (B, "<", $vcf_B) || die " ERROR: Unable to open $vcf_B\n";
    while (my $line=<B>) {
        chomp $line;
        next if $line =~/^#/;  
        my @tmp = split (/\t/, $line);
        my @info = split (";", $tmp[7]);
        next if $tmp[0] =~/_/;
        my $p = $tmp[1]+1;
	    my $var = "$tmp[3]/$tmp[4]";
	    #$vcf{"$tmp[0]\t$tmp[1]\t$p"}{B} = $Genotypes{$var};
        $vcf{"$tmp[0]\t$tmp[1]\t$p"}{B} = 1;		
    }
    close B;
    my $counter = 0;
    open (REGIONS, "<", $regions) || die " ERROR: Unable to open $regions\n";
    while (my $line=<REGIONS>) {
        chomp $line;
        my @tmp   = split (/\t/, $line);
        my $chr   = $tmp[0];
        my $start = $tmp[1];
        my $end   = $tmp[2];
        for (my $i = $start; $i <= $end; $i++) {
	        $counter++;
	        last if $counter > 200000;
            my $p = $i+1;
            foreach my $sample ( sort keys %samples ) {
                my $position = "$chr\t$i\t$p";
                if (!exists $vcf{"$chr\t$i\t$p"}{A}) {
                    push @array_A, "0";
                }
                else {
                    my $position = "$chr\t$i\t$p";
                    push @array_A, $vcf{$position}{A};
                }	
                if (!exists$vcf{"$chr\t$i\t$p"}{B}) {
                    push @array_B, "0";
                }
                else {
                    push @array_B, $vcf{$position}{B};
                }
            }
        }
    }
    #print " Calculating correlation\n";
    my $cor = correlation( \@array_A, \@array_B );
    $cor=~s/,/./; 
    return $cor;
}

##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -v,--vcf_dir   STRING   Input directory with VCF files
 -b,--bed_dir   STRING   Input directory where with BED files ending with *shared.bed
 -o,--outdir    STRING   Output directory
 -t,--threads   INT      Number of CPU cores (default = 4)\n\n";
	exit;
}