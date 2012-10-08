#
# AssemblyLine: transcriptome meta-assembly from RNA-Seq
#
# Copyright (C) 2012 Matthew Iyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# load libraries
library(rpart)
library(ROCR)

# parse command line
args <- commandArgs(trailingOnly=TRUE)
inputFile <- args[1]
outputPrefix <- args[2]

# output files
resultFile <- paste(outputPrefix, ".classify.txt", sep="")
cutoffFile <- paste(outputPrefix, ".cutoffs.txt", sep="")
treePlotFile <- paste(outputPrefix, ".ctree.pdf", sep="")
cpPlotFile <- paste(outputPrefix, ".cpplot.pdf", sep="")
prPlotFile <- paste(outputPrefix, ".prcurve.pdf", sep="");
cvfPlotFile <- paste(outputPrefix, ".cutoff_vs_f.pdf", sep="")
rocPlotFile <- paste(outputPrefix, ".roc.pdf", sep="")
aucFile <- paste(outputPrefix, ".auc.txt", sep="")

# read input data
t <- read.table(inputFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# setup classification tree parameters
minsplit <- max(30, 0.001 * nrow(t))
minbucket <- minsplit / 3.0
maxdepth <- 20
cp <- 0.0001
xval <- 10
control <- rpart.control(minsplit=minsplit,minbucket=minbucket,maxdepth=maxdepth,cp=cp,xval=xval,usesurrogate=2)

# setup case weights to balance ratio of observations in each class
total_obs <- nrow(t)
case_weights <- rep(1.0, nrow(t))
x <- table(t$annotated)
if (x["0"] > 0) {
	case_weights[which(t$annotated == 0)] <- (x["1"] / x["0"])
}

# run classification
fit <- rpart(annotated ~ length + num_exons + score + mean_score + mean_recurrence, 
		data=t, weights=case_weights, method="class", control=control)

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
pred <- predict(fit.prune, data=t, type="prob")
if (!("1" %in% colnames(pred))) {
	# there are no "1" in predicted probs so just
	# set all results to zero
	pred <- rep(0.0, dim(t)[1])
} else {
	pred <- pred[,"1"]
}

# compute statistics at different cutoff values
cutoffs <- seq(0, 1, by=0.001)
cols <- c("cutoff", "tp", "fn", "fp", "tn", "prec", "rec", "spec", "F", "balacc")
cutoffTable <- data.frame(matrix(nrow=length(cutoffs),ncol=length(cols), dimnames=list(NULL, cols)))
for (i in 1:length(cutoffs)) {
	cutoff <- cutoffs[i]
	tp <- length(which((pred >= cutoff) & (t$annotated == 1)))
	fn <- length(which((pred < cutoff) & (t$annotated == 1)))
	fp <- length(which((pred >= cutoff) & (t$annotated == 0)))
	tn <- length(which((pred < cutoff) & (t$annotated == 0)))
	prec <- 0.0
	rec <- 0.0
	spec <- 0.0
	F <- 0.0
	if ((tp + fp) > 0) {
		prec <- tp / (tp + fp)
	} 
	if ((tp + fn) > 0) {
		rec <- tp / (tp + fn)
	}
	if ((tn + fp) > 0) {
		spec <- tn / (tn + fp)
	}
	if ((prec + rec) > 0) {		
		F <- 2 * (prec * rec) / (prec + rec)
	}
	balacc <- (rec + spec) / 2.0
	cutoffTable[i,] <- c(cutoff, tp, fn, fp, tn, prec, rec, spec, F, balacc)
}

# write prediction output file
resultTable <- cbind(t, pred)
write.table(resultTable, file=resultFile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
# write statistics output file
write.table(cutoffTable, file=cutoffFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# print statistics from classification
if (fit.ok) {
	print(fit)
	printcp(fit)
	print(fit.prune)
	printcp(fit.prune)
}

# create plots
# attractive plot of tree
pdf(file=treePlotFile)	
try(plot(fit.prune, uniform=TRUE))
try(text(fit.prune, use.n=TRUE, all=TRUE, cex=0.7))
dev.off()
# CP plot
pdf(file=cpPlotFile)	
try(plotcp(fit))
dev.off()

# ROCR plots
pred.obj <- prediction(pred, t$annotated);
# precision-recall curve         
pdf(file=prPlotFile)
perf.rp <- performance(pred.obj, "prec", "rec");
try(plot(perf.rp, colorize=TRUE, main="Precision vs. Recall"));
dev.off()
# cutoff vs F measure
pdf(file=cvfPlotFile) 
perf.f <- performance(pred.obj, "f")
try(plot(perf.f, main="Cutoff vs. F-measure"));
dev.off()
# ROC
pdf(rocPlotFile)
perf.roc <- performance(pred.obj, "tpr", "fpr")
try(plot(perf.roc, colorize=TRUE, main="ROC"))
dev.off()
# AUC
auc.tmp <- performance(pred.obj, "auc");
auc <- as.numeric(auc.tmp@y.values)
write.table(auc, file=aucFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
