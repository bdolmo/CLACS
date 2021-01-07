#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $infile = $ARGV[0];

if (@ARGV < 2 ) {
	Help();
}
 my $devNull = ">/dev/null 2>&1";

 my $inputFile;
 my $outputDir;
 
 Help () if (@ARGV < 2 or !GetOptions(
	'input|i=s'=>\$inputFile,
	'outdir|o=s'=>\$outputDir
	)
 );

if (!-e $inputFile) {
   print " ERROR: $inputFile does not exist\n"; exit;
}

my $Rscript = `which Rscript`;
chomp $Rscript;
if (!$Rscript){
    print " ERROR: Rscript is not available on PATH\n";
    exit;
} 

createPlot($inputFile, $outputDir);

############
sub createPlot {
    my $inputFile = shift;
    my $outputDir= shift;
	my $libraries = "library(RColorBrewer)\nlibrary(gplots)\n";

	open (R, ">", "$outputDir/plotCorrelationMatrix.R") || die " ERROR: Cannot open plotCorrelationMatrix.R\n";
	print R "$libraries\n";
		
	print R "mydata<-read.csv(file=\"$inputFile\", sep =\'\t\', check.names = FALSE, header=TRUE)\n";
	print R "png(\"$outputDir/Heatmap.png\", res = 300, height=2800, width=2800)\n";
    print R "mat_data <- data.matrix(mydata[,2:ncol(mydata)])\n";
    print R "rnames <- mydata[,1]\n";
    print R "rownames(mat_data) <- rnames\n";
	#print R "my_palette <- colorRampPalette(c(\"white\",\"#327aff\",\"#6f6fff\", \"#196aff\", \"darkblue\"))(n = 150)\n";
	print R "my_palette <- brewer.pal(n = 9, name = \"Blues\")\n";

	# Creating correlation matrix
	print R "heatmap.2(mat_data, col=my_palette, sepcolor=\"black\", trace=\"none\",margins =c(12,9))\n";
	print R "dev.off()\n";
	close R;

    my $cmd = "$Rscript $outputDir/plotCorrelationMatrix.R $devNull";
    system $cmd;
}

############
sub Help {
	print " \n Usage: $0 <options>
 Options:
 -i,--input    STRING   Input TSV file
 -o,--outdir   STRING   Output directory\n\n";
	exit;
}


