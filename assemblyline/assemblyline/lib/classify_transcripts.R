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
fit <- rpart(sense_overlap ~ length + num_exons + intronic_overlap + mydensity + recur + avgdensity, data=t, method="anova", control=control)

# find optimal cp to prune
prune_cp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]

# prune tree
fit.prune <- prune(fit, cp=prune_cp)

# get predicted classes
result <- predict(fit.prune, data=t, type="vector")

# print statistics from classification
print(fit.prune)
printcp(fit)
printcp(fit.prune)

# create postscript plots
postscript(file=epsfile)

# create "attractive" plot of tree
post(fit.prune, title="Transcript Classification Tree")

# save histogram of distribution of the predictions
hist(result, breaks=50)

# determine an optimal cutoff value
sense_overlap_threshold=0.75
cutoffs <- seq(0, 1, by=0.001)
acc <- NA
i <- 1
#print("Performance at different cutoffs")
#print("================================")
for (cutoff in cutoffs) {
	tp <- length(which((result >= cutoff) & (t$sense_overlap >= sense_overlap_threshold)))
	fn <- length(which((result < cutoff) & (t$sense_overlap >= sense_overlap_threshold)))
	fp <- length(which((result >= cutoff) & (t$sense_overlap < sense_overlap_threshold)))
	tn <- length(which((result < cutoff) & (t$sense_overlap < sense_overlap_threshold)))
	spec <- tn / (tn + fp)
	sens <- tp / (tp + fn)
	F <- 2 * (spec * sens) / (sens + spec)
	acc[i] <- F
	i <- i + 1
	#pvals[i] <- fisher.test(matrix(c(tp,fn,fp,tn), nrow=2))$p.value
	#print(paste("CUTOFF", cutoff, "TP", tp, "FN", fn, "FP", fp, "TN", tn, "NEGLOGP", -log10(pvals[i])))
}
plot(cutoffs, acc)
dev.off()
best_cutoff <- cutoffs[min(which(acc == max(acc)))]
write(best_cutoff, file=cutoff_file, ncolumns=1) 

# write to output file
m <- matrix(result)
rownames(m) <- t$tx_id
write.table(m, file=outfile, quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")
