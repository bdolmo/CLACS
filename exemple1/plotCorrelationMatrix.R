library(RColorBrewer)
library(gplots)

mydata<-read.csv(file="exemple1/sample_correlations.tsv", sep ='	', check.names = FALSE, header=TRUE)
png("exemple1//Heatmap.png", res = 300, height=2800, width=2800)
mat_data <- data.matrix(mydata[,2:ncol(mydata)])
rnames <- mydata[,1]
rownames(mat_data) <- rnames
my_palette <- brewer.pal(n = 9, name = "Blues")
heatmap.2(mat_data, col=my_palette, sepcolor="black", trace="none",margins =c(12,9))
dev.off()
