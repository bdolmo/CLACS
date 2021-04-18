#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Config;

my $dirname = dirname(__FILE__);

my $bwa;
my $star;
my $samtools;
my $devNull = ">/dev/null 2>&1";

# Setting binaries for either mac (darwin) or Linux systems
if ($Config{osname} =~/darwin/) {
	$bwa = ( -x "$dirname/bin/darwin/bwa" )
	? "$dirname/bin/darwin/bwa"
	:  die " ERROR: Unable to execute bwa";

	$star = ( -x "$dirname/bin/darwin/STAR" )
	? "$dirname/bin/darwin/STAR"
	:  die " ERROR: Unable to execute STAR";

	$samtools = ( -x "$dirname/bin/darwin/samtools" )
	? "$dirname/bin/darwin/samtools"
	:  die " ERROR: Unable to execute samtools";
}
else {
	$bwa = ( -x "$dirname/bin/linux/bwa" )
	? "$dirname/bin/linux/bwa"
	:  die " ERROR: Unable to execute bwa";

	$star = ( -x "$dirname/bin/linux/STAR" )
	? "$dirname/bin/linux/STAR"
	:  die " ERROR: Unable to execute STAR";

	$samtools = ( -x "$dirname/bin/linux/samtools" )
	? "$dirname/bin/linux/samtools"
	:  die " ERROR: Unable to execute samtools";
}

# Check if java is available
my $java = `which java`; chomp $java;
if (!$java) {
	print " ERROR: java was not found on PATH\n";
	exit;
}
my $picard = "$dirname/bin/common/picard.jar";

my $inputDir;
my $outputDir;
my $reference;
my $mode = "dna";
# Using 4 threads in default 
my $threads = 4;
my $genomeVersion = "hg19";

Help () if (@ARGV < 1 or !GetOptions(
	'input|i=s'=>\$inputDir,
	'outdir|o=s'=>\$outputDir,
	'mode|m=s'=>\$mode, 
	'genome|g=s'=>\$reference,
	'reference|r=s'=>\$genomeVersion,
	'threads|t=i'=>\$threads
	)
);

if (! $inputDir || !-e $inputDir) {
    print " ERROR: Please, introduce a valid input directory\n";
	exit;
}

if ($mode !~/rna/i && $mode !~/dna/){
	print " ERROR: --mode parameter does only accept: rna or dna\n";
	exit;
} 

if (!-e $reference){
	print " ERROR: Missing a reference genome in FASTA format\n";
	exit;
}
if (!$genomeVersion) {
	print " ERROR: Missing genome version (hg19 or hg38)\n";
	Help();
}
if ($genomeVersion ne 'hg19' && $genomeVersion ne 'hg38') {
	print " ERROR: Incorrect genome version. Please, introduce hg19 or hg38\n";
	Help();
}

my $refDir = dirname($reference);

# Create output directory
mkdir $outputDir;

# Gathering paired fastq files and mapping
my @FASTQS = glob ("$inputDir/*.fastq.gz");

# Check if alternative fq suffix
if (!@FASTQS) {
	@FASTQS = glob ("$inputDir/*.fq.gz");
}
if (!@FASTQS) {
	print " ERROR: No gzipped FASTQ files were found on $inputDir input directory\n";
	exit;
}

# Checking all mapper indexes are present, if not, perform index generation 
if ($mode =~/dna/i) {
	# BWA
	createBwaIndex($reference);
}
elsif ($mode =~/rna/i) {
	# START
	createStartIndex($reference);
} 

my $localTime = localtime();
print "$localTime - INFO: bwa found: $bwa\n";
print "$localTime - INFO: samtools found: $samtools\n";
print "$localTime - INFO: picard found: $picard\n";

my %seenFastq = ();

my $mappingLog    = "$outputDir/Mapping.log";
my $duplicatesLog = "$outputDir/Duplicates.log";
open (LOG, ">", $mappingLog) || die " ERROR: Unable to open $mappingLog\n";
open (DUPLICATES, ">", $duplicatesLog) || die " ERROR: Unable to open $duplicatesLog\n";
print DUPLICATES "SAMPLE\tTOTAL_READS\tPCR_DUPLICATES\t\%DUPLICATES\n";

foreach my $fq (@FASTQS) {

	# Get the sample name 
    my $name = basename ( $fq );
    $name = ( split /_/, $name )[0] if $name =~/_/;

	$name =~s/.fastq.gz// if $name =~/.fastq.gz/;
	$name =~s/.fq.gz// if $name =~/.fq.gz/;
	$name =~s/.fastq// if $name =~/.fastq/;
	$name =~s/.fq// if $name =~/.fq/;

	$seenFastq{$name}++;
	next if ($seenFastq{$name} > 1);

    my $fastq1 = $fq;
    my $fastq2 = $fastq1;

	if ($fastq1 =~/R1/) {
		$fastq2 =~s/R1/R2/;
	}
	elsif ($fastq1 =~/R2/) {
		$fastq2 =~s/R2/R1/;
	}
	if ($fastq1 !~/R1/ or $fastq1 !~/R2/) {
		$fastq2 = "";
	}

	my $localTime = localtime();
	print "$localTime - INFO: mapping and sorting $name sample\n";
	print LOG "$localTime - INFO: mapping and sorting $name sample\n";

	my $cmd;
	# For paired-end
	if (-e $fastq1 && -e $fastq2 ) {
		if ($mode =~/dna/i ) { 
			$cmd = "$bwa mem -t $threads -R \'\@RG\\tID:$name\\tSM:$name\' $reference $fastq1 $fastq2";
			$cmd .= " | grep -v -e 'XA:Z:' -e 'SA:Z:' | $samtools sort -O BAM -o $outputDir/$name.bam - $devNull";
			system $cmd if !-e "$outputDir/$name.bam";
		} 
		if ($mode =~/rna/i ) { 
			$cmd = "$star --genomeDir $refDir --readFilesIn $fastq1 $fastq2 --readFilesCommand zcat --outSAMmapqUnique 60 --outFilterScoreMinOverLread 0.33  --outFilterMatchNminOverLread 0.33  --runThreadN $threads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outputDir/$name";
		} 
	}
	# For single-end
    if ( -e $fastq1 && !-e $fastq2 ) {
		if ($mode =~/dna/i ) { 
			$cmd = "$bwa mem -t $threads -R \'\@RG\\tID:$name\\tSM:$name\' $reference $fastq1";
			$cmd .= "| grep -v -e 'XA:Z:' -e 'SA:Z:' | $samtools sort -O BAM -o $outputDir/$name.bam - $devNull";
			system $cmd if !-e "$outputDir/$name.bam";
		}
		if ($mode =~/rna/i ) { 
			$cmd = "$star --genomeDir $refDir --readFilesIn $fastq1 --readFilesCommand zcat --outSAMmapqUnique 60 --outFilterScoreMinOverLread 0.33  --outFilterMatchNminOverLread 0.33  --runThreadN $threads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outputDir/$name";
		} 

	}
	$localTime = localtime();
	print "$localTime - INFO: removing PCR duplicates from sample $name\n";
	print LOG "$localTime - INFO: removing PCR duplicates from sample $name\n";

	if ( -e "$outputDir/$name.bam" ) {

		# Index raw bam file
		my $cmd = "$samtools index $outputDir/$name.bam";
		system $cmd;

		# Remove duplicates with picard
		$cmd = "$java -jar $picard MarkDuplicates INPUT=$outputDir/$name.bam OUTPUT=$outputDir/$name.nodup.bam METRICS_FILE=$outputDir/$name.DUPLICATES.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT $devNull";
		system $cmd if !-e "$outputDir/$name.nodup.bam";

		unlink("$outputDir/$name.bam");
		unlink("$outputDir/$name.bam.bai");

		# Index bam file without duplicates
		$cmd = "$samtools index $outputDir/$name.nodup.bam";
		system $cmd;

		rename "$outputDir/$name.nodup.bam", "$outputDir/$name.bam";
		rename "$outputDir/$name.nodup.bam.bai", "$outputDir/$name.bam.bai";
	}
	else {
		print "$localTime - ERROR: non-existent $outputDir/$name.bam\n";
		next;
	}
	$localTime = localtime();
	print "$localTime - INFO: counting total reads of sample $name\n";
	print LOG "$localTime - INFO: counting total reads of sample $name\n";

	my $total_reads = `$samtools view -F 0x4 $outputDir/$name.nodup.bam | cut -f 1 | sort | uniq | wc -l`;
	chomp $total_reads;

	my $perc_duplicates = `cat $outputDir/$name.DUPLICATES.txt | grep 'PERCENT_DUPLICATION | tail -3 | head -1 | cut -f 9`;
	chomp $perc_duplicates;

	print DUPLICATES "$name\t$total_reads\t$perc_duplicates\n";
}
close LOG;
close DUPLICATES;

##################################
sub createStarIndex {
	my $reference = shift;

	# Check that all bwa index files are available
	my $refDir = dirname($reference);
	my $localTime = localtime();
	my $flag = 0;

	my @indexFiles = ( 
		"chrLength.txt",
		"chrNameLength.txt",
		"chrName.txt",
		"chrStart.txt",
		"exonGeTrInfo.tab",
		"exonInfo.tab",
		"geneInfo.tab",
		"Genome",
		"genomeParameters.txt",
		"SA",
		"SAindex",
		"sjdbInfo.txt",
		"sjdbList.fromGTF.out.tab",
		"sjdbList.out.tab",
		"transcriptInfo.tab"
	);
	foreach my $file (@indexFiles) {
		if (!-e "$refDir/$file") {
			$flag = 1;
		}
	}

	if ($flag == 1) {

		my $doIndex = 0;

		while (1) {
			print "\nIndex files for STAR were not detected\n";
			print "Do you wish to create a new STAR index?: If yes enter \'y\' if no enter \'n\'\n";
			my $response = <STDIN>;
			chomp $response;
			if ($response eq 'y') {
				$doIndex = 1;
				last;
			}
			elsif ($response eq 'n') {
				$doIndex = 0;
				exit;
				last;
			}
			else {
				print " Incorrect option $response. Please, enter \'y\' or \'n\'\n";
			}
		}
		if ($doIndex) {
			# Download GTF files for human transcripts
			my $gencodeGtf = "gencode.v37.annotation.gtf.gz";
			if ($genomeVersion eq 'hg19') {
				my $cmd = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/$gencodeGtf?dl=0 -O $refDir/$gencodeGtf";
				system $cmd;
			}
			else {
				$gencodeGtf = "gencode.v37lift37.annotation.gtf.gz";
				my $cmd = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/$gencodeGtf?dl=0 -O $refDir/$gencodeGtf";
				system $cmd;
			} 
			my $cmd = "$star --runThreadN $threads --runMode genomeGenerate --genomeDir $refDir --genomeFastaFiles $reference --sjdbGTFfile $refDir/$gencodeGtf";
			system $cmd;
		}
	}
	if ($flag == 0) {
		print "$localTime - INFO: All STAR index files are available\n";
	}
} 


##################################
sub createBwaIndex {
	my $reference = shift;

	# Check that all bwa index files are available
	my $refDir = dirname($reference);
	my $localTime = localtime();

	my @indexFiles = (".amb", ".ann", ".pac", ".bwt", ".sa");
	my $flag = 0;
	foreach my $file (@indexFiles) {
		if (!-e "$reference$file") {
			$flag = 1;
		}
	}

	if ($flag == 1) {

		my $doIndex = 0;

		while (1) {
			print "\nIndex files for BWA were not detected\n";
			print "Do you wish to create a new BWA index?: If yes enter \'y\' if no enter \'n\'\n";
			my $response = <STDIN>;
			chomp $response;
			if ($response eq 'y') {
				$doIndex = 1;
				last;
			}
			elsif ($response eq 'n') {
				$doIndex = 0;
				exit;
				last;
			}
			else {
				print " Incorrect option $response. Please, enter \'y\' or \'n\'\n";
			}
		}
		if ($doIndex) {
			my $cmd = "$bwa index $reference";
			system $cmd;
		}
	}
	if ($flag == 0) {
		print "$localTime - INFO: All bwa index files are available\n";
	}

}

##################################
sub Help {
	print "\n Usage: $0 <options>
 Options:
 -i,--input      STRING   Input directory with gzipped FASTQ files
 -o,--outdir     STRING   Output directory
 -g,--genome     STRING   Reference genome in FASTA format
 -r,--reference  STRING   Genome version. Choose between [hg19, hg38] (default = hg19)
 -m,--mode       STRING   Alignment mode. Choose between [dna, rna] (default = dna)
 -t,--threads    INT      Number of CPU cores (default = 4)\n\n";
	exit;
}