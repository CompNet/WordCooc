#######################################################
# Additional processing for the tweet collection:
# - Read previously computed adjacency matrices.
# - Process inter-matrix distances.
# - Record distances.
# The process can be performed in a parallel way.
# 
# Vincent Labatut 3/2015
#
# setwd("D:/Eclipse/workspaces/Extraction/")
# setwd("~/Eclipse/workspaces/Extraction/")
# source("WordCooc/src/tweets/net-distance.R")
#######################################################
library("parallel")		# parallel packages
library("foreach")
library("doParallel")

source("WordCooc/src/common/misc.R")


# Whether the co-occurrence graphs should be directed (TRUE)
# or not (FALSE)
directed <- FALSE

# set up in/out folders
#in.folder <- "WordCooc/in/tweets/"
#out.folder <- "WordCooc/out/tweets/"
#log.file <- "WordCooc/out/log.txt"
##
#in.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/training/out_vincent/"
#out.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/training/out_vincent_features/"
##
#in.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/test/out_vincent/"
#out.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/test/out_vincent_features/"
##
#in.folder <- "~/work/data/training/"
#out.folder <- "~/work/data/training_features/"
#log.file <- "~/work/data/training_features/log.txt"
##
in.folder <- "~/work/data/test/"
out.folder <- "~/work/data/test_features/"
log.file <- "~/work/data/test_features/log.txt"

# set up log
dir.create(out.folder,recursive=TRUE,showWarnings=FALSE)
sink(file=log.file, append=FALSE, split=TRUE)

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

# set up parallel processing
core.nb <- detectCores()
par.clust <- makeCluster(core.nb/2)
registerDoParallel(par.clust)

##################################################################
# process inter-user distances based on the graphs
cat("Start processing the distances\n")
# treat each pair of users (i.e. graphs)
#for(i in 1:(length(text.files)-1))
foreach(i=1:(length(text.files)-1)) %dopar%
{	# set up folder
	text.file1 <- text.files[i]
	name1 <- strsplit(text.file1,"_texts.tsv.out")[[1]]
	subfolder1 <- paste(out.folder,text.file1,"/",sep="")
	cat("..Processing file ",text.file1," (",i,"/",(length(text.files)-1),")\n",sep="")
	
	# init result vectors
	result <- c()
	names1 <- c()
	names2 <- c()
	
	# get (partial) adjacency matrix
#	m <-
	co.count1 <- as.matrix(read.table(file=paste(subfolder1,"coocurrences.txt",sep=""),header=TRUE,row.names=1,
					check.names=FALSE)) # otherwise, an "X" is added in front of the colnames
#	co.count1 <- matrix(0,nrow=length(terms),ncol=length(terms))
#	rownames(co.count1) <- terms
#	colnames(co.count1) <- terms
#	co.count1[rownames(m),colnames(m)] <- m
	
	# compare to all other matrices (located after this one in the list of files)
	for(j in (i+1):length(text.files))
	{	# retrieve the other matrix
		text.file2 <- text.files[j]
		name2 <- strsplit(text.file2,"_texts.tsv.out")[[1]]
		subfolder2 <- paste(out.folder,text.file2,"/",sep="")
		cat("....Versus file ",text.file2," (",(j-i),"/",(length(text.files)-i),")\n",sep="")
		
#		m <- 
		co.count2 <- as.matrix(read.table(file=paste(subfolder2,"coocurrences.txt",sep=""),header=TRUE,row.names=1,check.names=FALSE))
#		co.count2 <- matrix(0,nrow=length(terms),ncol=length(terms))
#		rownames(co.count2) <- terms
#		colnames(co.count2) <- terms
#		co.count2[rownames(m),colnames(m)] <- m
		
		common.terms <- sort(unique(c(rownames(co.count1),rownames(co.count2))))
		# first matrix
		m1 <- matrix(0,nrow=length(common.terms),ncol=length(common.terms))
		rownames(m1) <- common.terms
		colnames(m1) <- common.terms
		m1[rownames(co.count1),colnames(co.count1)] <- co.count1
		# second matrix
		m2 <- matrix(0,nrow=length(common.terms),ncol=length(common.terms))
		rownames(m2) <- common.terms
		colnames(m2) <- common.terms
		m2[rownames(co.count2),colnames(co.count2)] <- co.count2
		
		# process distance and add to result vectors
		# taken from http://math.stackexchange.com/questions/507742/distance-similarity-between-two-matrices
		m <- m1 - m2
		d <- sqrt(sum(m*m))
		# note: for undirected networks, the lower triangle should be ignored when summing
		# however, here the trace is zero, so this doesn't affect the final result
		result <- c(result,d)
		names1 <- c(names1,name1)
		names2 <- c(names2,name2)
		
		# record vectors
		df <- data.frame("User1"=names1, "User2"=names2, "Distance"=result)
		write.table(x=df,file=paste(subfolder1,"distances.txt",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE)
	}
}	
## record vectors
#df <- data.frame("User1"=names1, "User2"=names2, "Distance"=result)
#write.table(x=df,file=paste(out.folder,"distances.txt",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE)

# stop parallel mode
stopCluster(par.clust)

# disable logging
sink()
