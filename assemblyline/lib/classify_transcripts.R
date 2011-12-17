# parse command line
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
cutoff_file <- args[3]
plotfile <- args[4]

# load libraries
library(rpart)

# read input data
t <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# setup classification tree parameters
minsplit <- max(30, 0.001 * nrow(t))
minbucket <- minsplit / 3.0
maxdepth <- 20
cp <- 0.0001
xval <- 10
control <- rpart.control(minsplit=minsplit,minbucket=minbucket,maxdepth=maxdepth,cp=cp,xval=xval,usesurrogate=2)

# run classification
fit <- rpart(annotated ~ length + num_exons + cov + fpkm + recur + avgdensity, data=t, method="class", control=control)

# check whether classification worked
fit.ok <- !is.nan(fit$cptable[1,1])

if (dim(fit$cptable)[1] == 1) {
	fit.prune <- fit
} else {
	# find optimal cp to prune
	prune_cp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
	# prune tree
	fit.prune <- prune(fit, cp=prune_cp)
}

# get predicted probabilities
result <- predict(fit.prune, data=t, type="prob")
if (!("1" %in% colnames(result))) {
	# there are no "1" in predicted probs so just
	# set all results to zero
	result <- rep(0.0, dim(t)[1])
} else {
	result <- result[,"1"]
}

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

# write to output file
m <- matrix(cbind(t$chrom, t$start, t$tx_id, t$annotated, t$recur, t$avgdensity, result), ncol=7)
write.table(m, file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# print statistics from classification
if (fit.ok) {
	print(fit)
	printcp(fit)
	print(fit.prune)
	printcp(fit.prune)
}

# create plots
pdf(file=plotfile)	
# create attractive plot of tree
try(plot(fit.prune, uniform=TRUE))
try(text(fit.prune, use.n=TRUE, all=TRUE, cex=0.7))
# show CP plot
try(plotcp(fit))
# plot histogram of distribution of the probabilities
try(hist(result, breaks=50))	
# plot overall cutoff
try(plot(cutoffs, F, main="Cutoff vs. F-measure"))	
# TODO: ROC curve
dev.off()	
