#
# AssemblyLine: transcriptome meta-assembly from RNA-Seq
#
# Copyright (C) 2012,2013 Matthew Iyer
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
has_rocr <- require(ROCR, quietly=TRUE)
has_ggplot2 <- require(ggplot2, quietly=TRUE)

# define constants
SATURATION <- 1e-10

# define category types
SAME_STRAND = 0    
OPP_STRAND = 1
INTRONIC_SAME_STRAND = 2
INTRONIC_OPP_STRAND = 3
INTRONIC_AMBIGUOUS = 4
INTERLEAVING = 5
INTERGENIC = 6
# two major category classes
INTRONIC_LIKE = c(INTRONIC_SAME_STRAND, INTRONIC_AMBIGUOUS)
INTERGENIC_LIKE = c(OPP_STRAND, INTRONIC_OPP_STRAND, INTERLEAVING, INTERGENIC)

bandwidth.nrd <- function (x) {
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2L] - r[1L])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}

kde2d <- function (x, y, h, n = 25, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) 
        stop("only finite values are allowed in 'lims'")
    n <- rep(n, length.out = 2L)
    gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
    gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
    h <- if (missing(h)) 
        c(bandwidth.nrd(x), bandwidth.nrd(y))
    else rep(h, length.out = 2L)
    if (any(h <= 0)) 
        stop("bandwidths must be strictly positive")
    h <- h/4
    ax <- outer(gx, x, "-")/h[1L]
    ay <- outer(gy, y, "-")/h[2L]
    z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), 
        , nx))/(nx * h[1L] * h[2L])
    list(x = gx, y = gy, z = z)
}

interp.surface <- function (obj, loc) {
    x <- obj$x
    y <- obj$y
    z <- obj$z
    nx <- length(x)
    ny <- length(y)
    lx <- approx(x, 1:nx, loc[, 1])$y
    ly <- approx(y, 1:ny, loc[, 2])$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    ex <- lx - lx1
    ey <- ly - ly1
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1
    return(z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + z[cbind(lx1 + 
        1, ly1)] * ex * (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 - 
        ex) * ey + z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}

kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=50,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,       # see option nlevels in contour
                      theta=30,         # see option theta in persp
                      phi=30,           # see option phi in persp
                      title="Title")
                      {
	z   <- d$z
	nrz <- nrow(z)
	ncz <- ncol(z)

	couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
	fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
	dim(fcol) <- c(nrz,ncz)
	fcol      <- fcol[-nrz,-ncz]

	par(mfrow=c(1,2),mar=c(0.5,0.5,0.5,0.5))
	persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,zlab="density")

	par(mar=c(2,2,2,2))
	image(d,col=couleurs)
	contour(d,add=T,nlevels=nlevels)
	box()
}

classify.kde2d <- function(ma, cl, clweights=c(1,1), kde2d.h=c(4,4), kde2d.n=100) {
	# find data ranges
	lims = c(range(ma[,1]), range(ma[,2]))
	# bivariate kernel density estimators
	has_cl0 <- sum(cl == 0) > 0
	if ( has_cl0 ) {
		x0 <- ma[cl==0,1]
		y0 <- ma[cl==0,2]
		d0 <- kde2d(x0,y0,n=kde2d.n,h=kde2d.h,lims=lims)
	} else {
		d0 <- NULL
	}
	x1 <- ma[cl==1,1]
	y1 <- ma[cl==1,2]
	d1 <- kde2d(x1,y1,n=kde2d.n,h=kde2d.h,lims=lims)
	x2 <- ma[cl==2,1]
	y2 <- ma[cl==2,2]
	d2 <- kde2d(x2,y2,n=kde2d.n,h=kde2d.h,lims=lims)
	# perform bilinear interpolation onto both surfaces
	z1 <- interp.surface(d1,ma)
	z2 <- interp.surface(d2,ma)
	# compute likelihood ratio at each observation
	lr <- (z2 + SATURATION)/(z1 + SATURATION)
	# recompute with weights
	zw1 <- z1 * clweights[1] + SATURATION
	zw2 <- z2 * clweights[2] + SATURATION
	lrw <- zw2 / zw1
	# return
	return(list(d0=d0,d1=d1,d2=d2,z=cbind(z1,z2),lr=lr,zw=cbind(zw1,zw2),lrw=lrw,
	       clweights=clweights))
}

classify.kde2d.plot <- function(fit, # object returned by classify.kde2d               
    numcolors=50,          # the number of colors to use
    zlim=c(min(z),max(z)), # limits in z coordinates
    numlevels=20,          # see option nlevels in contour
    theta=30,              # see option theta in persp
    phi=30)                # see option phi in persp
{
	# apply class weights
	z2 <- fit$clweights[2] * fit$d2$z + SATURATION
	z1 <- fit$clweights[1] * fit$d1$z + SATURATION
	z <- log10(z2 / z1)
	d <- list(x=fit$d2$x, y=fit$d2$y, z=z)

	# rest is from kde2dplot which modifications to allow
	# negative z values
	z   <- d$z
	nrz <- nrow(z)
	ncz <- ncol(z)

	colormap <- tail(topo.colors(trunc(1.4 * numcolors)),numcolors)
	zfrac <- (z - zlim[1])/(zlim[2]-zlim[1])
	fcol <- colormap[trunc(zfrac*(numcolors-1))+1]
	dim(fcol) <- c(nrz,ncz)
	fcol <- fcol[-nrz,-ncz]

	par(mfrow=c(1,2),mar=c(1,1,1,1))	
	persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,zlab="density")

	par(mar=c(2,2,2,2))
	image(d,col=colormap)
	contour(d,add=TRUE,nlevels=numlevels)
	box()
}

getPerformanceTable <- function(x, y) {
	pred.obj <- prediction(x, y)
	perf <- performance(pred.obj, "sens")
	cutoff <- perf@x.values[[1]]
	sens <- perf@y.values[[1]]
	spec <- performance(pred.obj, "spec")@y.values[[1]]
	#f <- performance(pred.obj, "f")@y.values[[1]]
	#acc <- performance(pred.obj, "acc")@y.values[[1]]
	balacc <- (sens + spec) / 2.0
	return(cbind(cutoff, sens, spec, balacc))
}

performanceAtCutoff <- function(cutoff, x, y) {
	tp <- length(which((x >= cutoff) & y))
	fn <- length(which((x < cutoff) & y))
	fp <- length(which((x >= cutoff) & (!y)))
	tn <- length(which((x < cutoff) & (!y)))
	sens <- 0.0
	spec <- 0.0
	f <- 0.0
	if ((tp + fn) > 0) {
		sens <- tp / (tp + fn)
	}
	if ((tn + fp) > 0) {
		spec <- tn / (tn + fp)
	}
	balacc <- (sens + spec) / 2.0
	return(c(tp, fp, fn, tn, sens, spec, balacc))
}

classifyAndWriteResults <- function(inp, 
	variables, 
	cl, 
	clweights, 
	prefix, 
	kde2d.h, 
	kde2d.n) 
{
	# classification
	fit <- classify.kde2d(inp[,variables], cl, clweights=clweights, kde2d.h=kde2d.h, kde2d.n=kde2d.n)
	log10lr <- log10(fit$lrw)	
	res <- data.frame(cbind(cl, log10lr))

	# performance of classification
	has_tests <- (sum(cl == 0) > 0)
	perfFile <- paste(prefix, ".perf.txt", sep="")
	cat("type", "auc", "cutoff", 
		"train.tp", "train.fp", "train.fn", "train.tn", "train.sens", "train.spec", "train.balacc",
		"test.tp", "test.fp", "test.fn", "test.tn", "test.sens", "test.spec", "test.balacc", "\n",
		file=perfFile, sep="\t")
	# training data
	res.train <- res[(cl == 1) | (cl == 2),]
	pred.obj <- prediction(res.train$log10lr, res.train$cl == 2)
	perf.train.roc <- performance(pred.obj, "tpr", "fpr")
	auc.train <- performance(pred.obj, measure = "auc")@y.values[[1]]
	perf.train.table <- getPerformanceTable(res.train$log10lr, res.train$cl == 2)	
	cutoff.train <- perf.train.table[which.max(perf.train.table[,"balacc"]),"cutoff"]
	cutoff.train.perf <- performanceAtCutoff(cutoff.train, res.train$log10lr, res.train$cl == 2)
	
	# test data
	if ( has_tests ) {
		res.test <- res[(cl == 0) | (cl == 1),]
		pred.obj <- prediction(res.test$log10lr, (res.test$cl == 0))
		perf.test.roc <- performance(pred.obj, "tpr", "fpr")
		auc.test <- performance(pred.obj, measure = "auc")@y.values[[1]]
		perf.test.table <- getPerformanceTable(res.test$log10lr, res.test$cl == 0)
		cutoff.test <- perf.test.table[which.max(perf.test.table[,"balacc"]),"cutoff"]
		cutoff.test.perf <- performanceAtCutoff(cutoff.test, res.test$log10lr, res.test$cl == 0)
		cutoff.train.ontest <- performanceAtCutoff(cutoff.train, res.test$log10lr, res.test$cl == 0)
		cutoff.test.ontrain <- performanceAtCutoff(cutoff.test, res.train$log10lr, res.train$cl == 2)
	} else {
		auc.test <- NA
		cutoff.test.perf <- rep(NA, 7)
		cutoff.train.ontest <- rep(NA, 7)
		cutoff.test.ontrain <- rep(NA, 7)
	}

	# performance
	cat("train", auc.train, cutoff.train, cutoff.train.perf, cutoff.train.ontest, "\n", file=perfFile, sep="\t", append=TRUE)
	cat("test", auc.test, cutoff.test, cutoff.test.ontrain, cutoff.test.perf, "\n", file=perfFile, sep="\t", append=TRUE)
	
	# result table
	pred.train <- (res[,"log10lr"] > cutoff.train)
	if ( has_tests ) {
		pred.test <- (res[,"log10lr"] > cutoff.test)
	} else {
		pred.test <- rep(NA, nrow(res))
	}
	res <- cbind(res, pred.train, pred.test)

	# plots
	pdfFile <- paste(prefix, ".plots.pdf", sep="")
	pdf(pdfFile)

	# density plots
	kde2dplot(fit$d1)
	kde2dplot(fit$d2)
	if ( has_tests ) {
		kde2dplot(fit$d0)
	}

	# log likelihood landscape plot
	classify.kde2d.plot(fit)
	# cutoff vs balanced accuracy
	par(mfrow=c(1,1),mar=c(5,5,2,2))
	plot(perf.train.table[,"cutoff"], perf.train.table[,"balacc"], pch=18, 
	     xlab="Cutoff", ylab="Balanced Accuracy", ylim=c(0.5,1), col="red")
	if ( has_tests ) {
		points(perf.test.table[,"cutoff"], perf.test.table[,"balacc"], pch=18, col="blue")
		legend("topright", c("train", "test"), fill=c("red","blue"))
	}

	# roc curves
	par(mfrow=c(1,1),mar=c(5,5,2,2))
	plot(perf.train.roc, col="red", main="ROC")
	if ( has_tests ) {
		plot(perf.test.roc, col="blue", add=TRUE)
		legend("bottomright", c("train", "test"), fill=c("red","blue"))
	}

	# density histogram
	p <- ggplot(res, aes(x=log10lr, fill=factor(cl)))
	p <- p + geom_density(alpha=0.2)
	p <- p + geom_vline(xintercept = cutoff.train, colour="red", linetype="longdash")
	if ( has_tests ) {
		p <- p + geom_vline(xintercept = cutoff.test, colour="blue", linetype="longdash")
	}
	plot(p)
	dev.off()

	return(res)
}

# parse command line
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
input_file <- paste(prefix, ".inp.txt", sep="")
min_obs <- 50
max_frac_test_obs <- 0.05
kde2d.n <- 50
kde2d.h <- c(5,5)
prior_mrna <- 0.85
prior_intronic <- 0.1
prior_intergenic <- 1 - (prior_mrna + prior_intronic)
variables <- c("mean_recurrence", "pctrank")

# read input data
tbl <- read.table(input_file, header=TRUE, sep="\t")

# setup test cases
num_mrna <- sum((tbl$category == SAME_STRAND) | (tbl$test == 1))
# limit the size of the test data to a fraction of the total data
frac_test <- sum(tbl$test == 1) / num_mrna
pct_test_to_use <- min(max_frac_test_obs, frac_test) / frac_test
num_tests <- round(pct_test_to_use * sum(tbl$test == 1))
whichtests <- sample(which(tbl$test == 1), num_tests)
# rewrite table with test case information
whichtests.categories <- tbl[whichtests, "category"]
tbl[tbl$test == 1, "category"] <- 0
tbl$test <- 0
tbl[whichtests, "test"] <- 1
tbl[whichtests, "category"] <- whichtests.categories

# divide known transcripts into training/test sets
train <- (tbl$category == 0) & (tbl$test == 0)
test <- (tbl$test == 1)
# divide unknown and test transcripts into classes
testknown <- (tbl$test == 1) & (tbl$category == SAME_STRAND)
testintronic <- (tbl$test == 1) & (tbl$category %in% INTRONIC_LIKE)
testintergenic <- (tbl$test == 1) & (tbl$category %in% INTERGENIC_LIKE)
intronic <- (tbl$test == 0) & (tbl$category %in% INTRONIC_LIKE)
intergenic <- (tbl$test == 0) & (tbl$category %in% INTERGENIC_LIKE)

# compute ratio of each type of transcript
frac_mrna <- num_mrna / nrow(tbl)
frac_intronic <- sum(intronic) / nrow(tbl)
frac_intergenic <- sum(intergenic) / nrow(tbl)

# compute weights based on priors
known_vs_intronic_ratio <- (frac_mrna / prior_mrna) / (frac_intronic / prior_intronic)
known_vs_intergenic_ratio <- (frac_mrna / prior_mrna) / (frac_intergenic / prior_intergenic)

# determine whether there are enough samples to perform classification
do_intronic <- (num_mrna >= min_obs) & (sum(intronic) >= min_obs)
do_intergenic <- (num_mrna >= min_obs) & (sum(intergenic) >= min_obs)

# write stats information
statsFile <- paste(prefix, ".info.txt", sep="")
cat("mrna", num_mrna, frac_mrna, "\n", file=statsFile, sep="\t")
cat("tests", num_tests, "\n", file=statsFile, append=TRUE, sep="\t")
cat("intronic", sum(intronic), frac_intronic, "\n", sep="\t", file=statsFile, append=TRUE)
cat("intergenic", sum(intergenic), frac_intergenic, "\n", sep="\t", file=statsFile, append=TRUE)
cat("prior_mrna_vs_intronic", known_vs_intronic_ratio, "\n", sep="\t", file=statsFile, append=TRUE)
cat("prior_mrna_vs_intergenic", known_vs_intergenic_ratio, "\n", sep="\t", file=statsFile, append=TRUE)
cat("do_intronic", do_intronic, "\n", sep="\t", file=statsFile, append=TRUE)
cat("do_intergenic", do_intergenic, "\n", sep="\t", file=statsFile, append=TRUE)

# classify intronic-like transcripts
log10lr.intronic <- rep(NA, nrow(tbl))
pred.train.intronic <- rep(NA, nrow(tbl))
pred.test.intronic <- rep(NA, nrow(tbl))
classes <- rep(NA, nrow(tbl))
if ( do_intronic ) {
	intronicrows <- (train | testintronic | intronic)
	inp <- tbl[intronicrows,]
	cl <- ifelse(testintronic, 0, ifelse(intronic, 1, 2))[intronicrows]
	clweights <- c(1, known_vs_intronic_ratio)
	intronicprefix <- paste(prefix, ".intronic", sep="")
	res <- classifyAndWriteResults(inp, variables, cl, clweights, intronicprefix, 
                                   kde2d.h=kde2d.h, kde2d.n=kde2d.n)
	log10lr.intronic[intronicrows] <- res[,"log10lr"]
	pred.train.intronic[intronicrows] <- res[,"pred.train"]
	pred.test.intronic[intronicrows] <- res[,"pred.test"]
}

# classify intergenic-like transcripts
log10lr.intergenic <- rep(NA, nrow(tbl))
pred.train.intergenic <- rep(NA, nrow(tbl))
pred.test.intergenic <- rep(NA, nrow(tbl))
if ( do_intergenic ) {
	intergenicrows <- (train | testintergenic | intergenic)
	inp <- tbl[intergenicrows,]
	cl <- ifelse(testintergenic, 0, ifelse(intergenic, 1, 2))[intergenicrows]
	clweights <- c(1, known_vs_intergenic_ratio)
	intergenicprefix <- paste(prefix, ".intergenic", sep="")
	res <- classifyAndWriteResults(inp, variables, cl, clweights, intergenicprefix, 
                                   kde2d.h=kde2d.h, kde2d.n=kde2d.n)
	log10lr.intergenic[intergenicrows] <- res[,"log10lr"]
	pred.train.intergenic[intergenicrows] <- res[,"pred.train"]
	pred.test.intergenic[intergenicrows] <- res[,"pred.test"]
}

# write results
tbl <- cbind(tbl, log10lr.intronic, pred.train.intronic, pred.test.intronic,
		     log10lr.intergenic, pred.train.intergenic, pred.test.intergenic)
resultsFile <- paste(prefix, ".out.txt", sep="")
write.table(tbl, file=resultsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
