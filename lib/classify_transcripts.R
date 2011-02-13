# parse command line
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
epsfile <- args[3]

# load libraries
library(rpart)

# read input data
t <- read.table(infile, header=TRUE, sep="\t")

# setup classification tree parameters
minsplit <- max(30, 0.001 * nrow(t))
minbucket <- minsplit / 3.0
maxdepth <- 20
cp <- 0.001
xval <- 10
control <- rpart.control(minsplit=minsplit,minbucket=minbucket,maxdepth=maxdepth,cp=cp,xval=xval)

# run classification
fit <- rpart(ann ~ length+num_exons+lanes+score, data=t, method="class", control=control)

# find optimal cp to prune
prune_cp <- fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]

# prune tree
fit.prune <- prune(fit, cp=prune_cp)

# print statistics from classification
print(fit)
printcp(fit)

# create attractive postscript plot of tree
post(fit, file=epsfile, title="Transcript Classification Tree")

# get predicted classes
result <- predict(fit, data=t, type="class")

# write to output file
write(as.vector(result), outfile, ncolumns=1)


