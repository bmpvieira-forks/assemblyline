# parse command line
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
cutoff_file <- args[3]
epsfile <- args[4]

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
fit <- rpart(known ~ length + num_exons + mydensity + recur + avgdensity, data=t, method="class", control=control)
print(fit)
printcp(fit)

if (dim(fit$cptable)[1] == 1) {
	fit.prune <- fit
} else {
	# find optimal cp to prune
	prune_cp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
	# prune tree
	fit.prune <- prune(fit, cp=prune_cp)
}

# print statistics from classification
print(fit.prune)
printcp(fit.prune)

# get predicted classes
result <- predict(fit.prune, data=t, type="prob")
result <- result[,"1"]

# create postscript plots
postscript(file=epsfile)

# create "attractive" plot of tree
#post(fit.prune, title="Transcript Classification Tree")

# save histogram of distribution of the predictions
hist(result, breaks=50)

# determine an optimal cutoff value
cutoffs <- seq(0, 1, by=0.001)
acc <- NA
i <- 1
for (cutoff in cutoffs) {
	# cutoff for all transcripts
	tp <- length(which((result >= cutoff) & (t$known == 1)))
	fn <- length(which((result < cutoff) & (t$known == 1)))
	fp <- length(which((result >= cutoff) & (t$known == 0)))
	tn <- length(which((result < cutoff) & (t$known == 0)))
	spec <- tn / (tn + fp)
	sens <- tp / (tp + fn)
	F <- 2 * (spec * sens) / (sens + spec)
	acc[i] <- F
	i <- i + 1
}

# plot overall cutoff
plot(cutoffs, acc, main="Overall accuracy")
dev.off()

best_cutoff <- cutoffs[min(which(acc == max(acc)))]
write(c(best_cutoff), file=cutoff_file, ncolumns=1) 

# write to output file
m <- matrix(cbind(t$known, t$within_gene, result), ncol=3)
rownames(m) <- t$tx_id
write.table(m, file=outfile, quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")
