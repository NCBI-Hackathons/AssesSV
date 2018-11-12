#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

upset_plotter = function(filename) {

	library(UpSetR)
	tp = read.csv(filename, header=T, sep=";")

	png("tp_concord_conflict.png", width = 1000, height = 500) # , width = 16, height = 8, paper = "USr", onefile = FALSE
	upset(tp, sets = colnames(tp)[2:ncol(tp)], order.by = "freq", sets.bar.color = "#56B4E9", point.size = 3, line.size = 1, text.scale=1.5)
	dev.off()

	fn = tp

	for (row in 1:nrow(fn)) {
		for (col in 2:ncol(fn)) {
			fn[row,col] = abs(1 - tp[row,col])
		}
		if (row %% 10000 == 0) 
			cat("Processed row", row, "\n")
	}

	png("fn_concord_conflict.png", width = 1000, height = 500) # , width = 16, height = 8, paper = "USr", onefile = FALSE
	upset(fn, sets = colnames(fn)[2:ncol(fn)], order.by = "freq", sets.bar.color = "#56B4E9", point.size = 3, line.size = 1, text.scale=1.5)
	dev.off()

}

# upset_plotter(args[1])