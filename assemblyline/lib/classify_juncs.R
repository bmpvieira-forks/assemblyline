#
# load libraries
#
library(rpart)

#
# define functions
#
classifyjuncs <- function(t, rpart_control, outfile, cutoff_file, plotfile) {		
	# run classification
	fit <- rpart(annotated ~ myreads + myoverhang + enough_reads + best_cov + best_overhang + motif + recur, data=t, method="class", control=rpart_control)
	
	# prune the tree
	if (dim(fit$cptable)[1] == 1) {
		fit.prune <- fit
	} else {
		# find optimal cp to prune
		prune_cp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
		# prune tree
		fit.prune <- prune(fit, cp=prune_cp)
	}
	
	# print statistics before/after pruning
	print("Tree Before pruning")
	print(fit)
	print("Tree After pruning")
	print(fit.prune)
	print("CP Before pruning")
	printcp(fit)
	print("CP After pruning")
	printcp(fit.prune)

	# get predicted probabilities
	#result <- predict(fit, data=t, type="class")
	result <- predict(fit.prune, data=t, type="prob")
	result <- result[,"1"]
	
	# determine an optimal cutoff value
	cutoffs <- seq(0, 1, by=0.001)
	F <- NA
	i <- 1
	for (cutoff in cutoffs) {
		tp <- length(which((result >= cutoff) & (t$annotated == 1)))
		fn <- length(which((result < cutoff) & (t$annotated == 1)))
		fp <- length(which((result >= cutoff) & (t$annotated == 0)))
		tn <- length(which((result < cutoff) & (t$annotated == 0)))
		if (tp == 0) {
			F[i] <- 0.0
		} else {
			prec <- tp / (tp + fp)
			rec <- tp / (tp + fn)
			F[i] <- 2 * (prec * rec) / (prec + rec)
		}
		i <- i + 1
	}
	
	# find cutoff with best F-measure of performance	
	best_cutoff <- cutoffs[min(which(F == max(F)))]
	write(c(best_cutoff), file=cutoff_file, ncolumns=1)
	
	# write to tab-delimited output file
	m <- cbind(t, result)
	write.table(m, file=outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
	
	# create plots
	pdf(file=plotfile)	
	# create attractive postscript plot of tree
	try(plot(fit.prune, uniform=TRUE))
	try(text(fit.prune, use.n=TRUE, all=TRUE, cex=0.7))
	# plot histogram of distribution of the probabilities
	try(hist(result, breaks=50))	
	# plot overall cutoff
	try(plot(cutoffs, F, main="Cutoff vs. F-measure"))	
	# TODO: ROC curve
	dev.off()	
}

# parse command line
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outprefix <- args[2]
outfile <- paste(outprefix, ".classify.juncs", sep="")
cutoff_file = paste(outprefix, ".classify.cutoff", sep="")
plotfile = paste(outprefix, ".classify.plots.pdf", sep="")

# read input data
dat <- read.table(infile, sep="\t", header=TRUE)

# setup classification tree parameters
cp <- 0.0001
minsplit <- max(30, cp * nrow(dat))
minbucket <- minsplit / 3.0
maxdepth <- 30
xval <- 10
control <- rpart.control(minsplit=minsplit,minbucket=minbucket,maxdepth=maxdepth,cp=cp,xval=xval,usesurrogate=2)

# run classification
classifyjuncs(dat, control, outfile, cutoff_file, plotfile)
