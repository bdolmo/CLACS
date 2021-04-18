#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Config;

my $dirname = dirname(__FILE__);

my $infile = $ARGV[0];

if (@ARGV < 2 ) {
	Help();
}
 my $devNull = ">/dev/null 2>&1";

 my $inputDir;
 my $outputDir;
 my $genomeVersion;
 my $regionsBed;
 my $normFactor = 10e7;

 Help () if (@ARGV < 3 or !GetOptions(
	'input|i=s'=>\$inputDir,
	'reference|r=s'=>\$genomeVersion,
    'regions|e=s'=>\$regionsBed,
	'normfactor|n=i'=>\$normFactor, 
	'outdir|o=s'=>\$outputDir
	)
 );

if (!-e $inputDir) {
   print " ERROR: $inputDir does not exist\n"; exit;
}
if (!-e $outputDir) {
    mkdir $outputDir;
}

if (!$genomeVersion) {
	print " ERROR: Missing genome version (hg19 or hg38)\n";
	Help();
}
if ($genomeVersion ne 'hg19' && $genomeVersion ne 'hg38') {
	print " ERROR: Incorrect genome version. Please, introduce hg19 or hg38\n";
	Help();
}

my $Rscript = `which Rscript`;
chomp $Rscript;
if (!$Rscript){
    print " ERROR: Rscript is not available on PATH\n";
    exit;
} 

my $bedtools;
my $samtools;
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

my @bamFiles = glob("$inputDir/*.bam"); 
if (!@bamFiles) {
	print " ERROR: No input bam files were found on $inputDir directory\n";
	exit;
}

my $bamStr =  join(" ", @bamFiles );
# Get read counts at specificed regions file
my $cmd = "$bedtools multicov -bams $bamStr -bed $regionsBed > $outputDir/raw.counts.bed";
system $cmd;

my $samples =  join("\t", map basename($_ ), @bamFiles);
$samples =~s/.bam//g; 

# Add header to file
$cmd = "sed -i \'1i chr\tstart\tend\tregion\t$samples\' $outputDir/raw.counts.bed";
system $cmd; 

my %HoS = (); 

for my $bam (@bamFiles){
    my $name = basename ($bam);
    $name =~s/.bam//;
	my $counts = `$samtools idxstats $bam | awk -F '\t\' \'{s+=\$3+\$4}END{print s}'`;
	chomp $counts;
    $HoS{$name}{COUNTS} = $counts;
    $HoS{$name}{SCALE}  = $normFactor/$counts;
} 

open (IN, "<", "$outputDir/raw.counts.bed" ) 
    || die " ERROR: Unable to open $outputDir/raw.counts.bed\n";  
open (OUT, ">", "$outputDir/normalized.counts.bed") 
    || die " ERROR: Unable to open $outputDir/normalized.counts.bed\n"; 

my @headerSamples = ();  
while (my $line=<IN>){
    chomp $line;
    my @tmp = split("\t", $line); 
    if ($line=~/^chr\tstart/){
        @headerSamples = @tmp[4..@tmp-1]  ;
        print OUT "$line\n";
        next;
    } 
    my $j = 0;
    my @normCountsArr = ();  
    for (my $i=4; $i<@tmp;$i++){
        my $sample = $headerSamples[$j];
        my $normCounts = $tmp[$i]*$HoS{$sample}{SCALE}; 
        push @normCountsArr, $normCounts;
        $j++;
    }
    print OUT join("\t", @tmp[0..3]) . "\t" . join("\t", @normCountsArr) . "\n";
} 

close IN;
close OUT;

plotHistogram("$outputDir/normalized.counts.bed", $outputDir);

##################################
sub plotHistogram {
    my $normCov   = shift;
    my $outputDir = shift;

    my @libs = (
        "library(RColorBrewer)",
        "library(ggplot2)",
        "library(tidyr)",
        "library(gtools)",
        "library(gridExtra)",
        "library(dplyr)",
        "library(grid)"
    );
    my $libraries = join("\n", @libs);

	open (R, ">", "$outputDir/plotHistogram.R") || die " ERROR: Cannot open $outputDir/plotHistogram.R\n";
	print R "$libraries\n";
	print R "mydata<-read.table(file=\"$normCov\", sep =\'\t\', check.names = FALSE, header=TRUE)\n";

    print R "mydata <- mydata %>% gather(Sample, Counts, 5:ncol(mydata))\n";
    print R "mydata\$region<-factor( mydata\$region, levels=unique(mixedsort(mydata\$region)))\n";
    print R "mydata\$region <- factor(mydata\$region,levels = c(\"A1\",\"A2\",\"A3\",\"A4\",\"A5\",\"A6\",\"A7\",\"A8\",\"A9\",\"A10\",\"A11\", \"A12\", \"A13\",\"GA1\", \"GA2\", \"GA3\", \"GB1\", \"GA4\", \"GB2\", \"GA5\", \"GB3\", \"GA6\", \"GA7\", \"GB4\", \"GA8\", \"GB5\", \"GA9\", \"GB6\", \"GA10\", \"GB7\", \"GA11\", \"GA12\"))\n";
    print R "pcdha <-c(\"A1\",\"A2\",\"A3\",\"A4\",\"A5\",\"A6\",\"A7\",\"A8\",\"A9\",\"A10\",\"A11\", \"A12\",\"A13\")\n";
    print R "pcdhb <-c(\"GA1\", \"GA2\", \"GA3\", \"GB1\", \"GA4\", \"GB2\", \"GA5\", \"GB3\", \"GA6\", \"GA7\", \"GB4\", \"GA8\", \"GB5\", \"GA9\", \"GB6\", \"GA10\", \"GB7\", \"GA11\", \"GA12\")\n";

    print R " j = 1\n";
    print R "
    p_list <- list()
    for (i in 1:14) {
        sample <- sample_list[i]
        out_png <- paste(\"$outputDir\", sample , \".png\", sep="")
        #png(out_png, res = 150, height=235, width=1600)
        sampledata<-mydata %>% filter(Sample == sample)
        data_pcdha<-filter(sampledata, region \%in\% pcdha)
        data_others<- filter(sampledata, region \%in\% pcdhb)

        min<- 0
        max1 <- max(data_pcdha\$Counts)
        max2 <- max(data_others\$Counts)
        max <- max1
        if (max2 > max) {
            max <- max2
        } 
        limits <- c(0, max)
        breaks <- seq(min, max, length.out=4)
        p_list[[j]] <- ggplot(data=data_pcdha, aes_string(x=data_pcdha\$region, y=data_pcdha\$Counts, fill=data_pcdha\$Counts)) + 
            geom_bar(stat=\"identity\",width=.75, colour=\"black\", size=0.25) + scale_fill_gradient(low=\"white\", high=\"lightseagreen\") +
            theme(plot.title=element_text(face=\"bold\")) + theme(legend.position=\"none\") +
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+ 
            theme(panel.background=element_rect(fill=\"white\", color=\"white\"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            theme(axis.line=element_blank()) + theme(plot.title=element_text(face=\"bold\")) + scale_y_continuous(limits=limits,breaks=breaks)
        p_list[[j+1]] <- ggplot(data=data_others, aes_string(x=data_others\$region, y=data_others\$Counts, fill=data_others\$Counts)) + 
            geom_bar(stat=\"identity\",width=.75, colour=\"black\", size=0.25) + scale_fill_gradient(low=\"white\", high=\"lightseagreen\") +
            theme(plot.title=element_text(face=\"bold\")) + theme(legend.position=\"none\") +
            theme(axis.text.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+ 
            theme(panel.background=element_rect(fill=\"white\", color=\"white\"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            theme(axis.line=element_blank()) + theme(plot.title=element_text(face=\"bold\")) + scale_y_continuous(limits=limits,breaks=breaks)
        p_list[[j]] <- ggplot_gtable(ggplot_build(p_list[[j]]))
        p_list[[j+1]] <- ggplot_gtable(ggplot_build(p_list[[j+1]]))
        p_list[[j]]\$heights <- p_list[[j+1]]\$heights
        j<-j+2
    }\n";
	print R "png(\"$outputDir/histogram.normalized.counts.png\", res = 250, height=3500, width=1600)\n";
    print R "grid.arrange(grobs=p_list, widths=c(0.4,0.56),ncol=2 )\n";
	print R "dev.off()\n";
	close R;

    my $cmd = "$Rscript $outputDir/plotHistogram.R $devNull";
    system $cmd;
} 

##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -i,--input       STRING   Input directory with *.normalized.bed.gz files
 -o,--outdir      STRING   Output directory
 -e,--regions     STRING   BED regions to extract counts
 -n,--normfactor  INT      Normalization factor (default=10e7)
 -r,--reference   STRING   Genome version. Choose between [hg19, hg38] (default = hg19)\n\n";
	exit;
}