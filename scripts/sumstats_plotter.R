#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

sumstats_plotter = function(sum_table, stat_type) {

	library(ggplot2)
	title = paste("Comparison of", stat_type)
	yaxis = paste("Number of", stat_type)
	filename = paste(stat_type,".png", sep = "")
	png(filename, width = 1000, height = 500)
	plot = ggplot(sum_table, aes(as.factor(sum_table[,1]), sum_table[,2])) + labs(title = title, y = yaxis, x = "Conditions") + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	print(plot)
	dev.off()

}

# summary_table = read.csv(args[1], header=T, sep="\t")
# sum_stat_types = c("False_Negatives", "False_Positives", "True_Positives", "F-score", "Precision", "Recall")
# for(sum_stat in 2:ncol(summary_table)) {
#	sumstats_plotter(summary_table[,c(1,sum_stat)], sum_stat_types[sum_stat - 1])
# }