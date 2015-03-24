#######################################################
# Main processing of the tweet collection:
# - Retrieve the text files.
# - Extract the co-occurrence networks.
# - Process a bunch of nodal measures.
# - Record the resulting networks.
# 
# One text file corresponds to all the tweets produced
# by a given user. Each line in the file represents a
# distinct tweet. This script generates one network
# for each user, containing the nodal topological measures.
# 
# Vincent Labatut 3/2015
#
# setwd("D:/Eclipse/workspaces/Extraction/")
# setwd("~/Eclipse/workspaces/Extraction/")
# source("WordCooc/extract-net-tweets.R")
#######################################################
library("igraph")
library("parallel")		# parallel packages
library("foreach")
library("doParallel")

source("WordCooc/com-measures.R")
source("WordCooc/misc.R")

# TODO for directed networks >> adapt the measures

# Whether or not to record secondary data such as co-occurrence
# networks, list of terms, co-occurrence matrix, etc.
record.secondary.data <- FALSE
# TRUE to record full vectors (including many zeros), i.e. a value
# for each word detected in the corpus, even if it does not appear
# anywhere in the considered user's tweets. FALSE allows focusing
# only on the words actually used.
output.full.matrix <- FALSE
# Whether the co-occurrence graphs should be directed (TRUE)
# or not (FALSE)
directed <- FALSE

# set up in/out folders
in.folder <- "WordCooc/in/tweets/"
out.folder <- "WordCooc/out/tweets/"
log.file <- "WordCooc/out/log.txt"
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
#in.folder <- "~/work/data/test/"
#out.folder <- "~/work/data/test_features/"
#log.file <- "~/work/data/test_features/log.txt"

# set up log
sink(file=log.file, append=FALSE, split=TRUE)

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

# set up parallel processing
core.nb <- detectCores()
par.clust <- makeCluster(core.nb/2)
registerDoParallel(par.clust)


# build the list of terms for the whole collection
terms.file <- paste(out.folder,"terms.txt",sep="")
if(file.exists(terms.file))
{	cat("Loading the corpus terms\n")
	terms <- as.matrix(read.table(terms.file))	
}else
{	cat("Identifying the corpus terms\n")
	terms <- c()
	for(i in 1:length(text.files))
	{	# read the file
		text.file <- text.files[i]
		cat("Reading sentences for file ",i,"/",length(text.files)," ",text.file,"\n",sep="")
		sentences <- readLines(paste(in.folder,text.file,sep=""))
		
		# remove empty sentences
		sentences <- sentences[nchar(sentences,type="bytes")>0]
		
		# split each line in columns using tab as a separator
		cat("Spliting sentences\n")
		cols <- strsplit(sentences, "\t")
		
		# split each sentence using space as a separator
		words <- lapply(cols,function(c) 
				{	sentence <- paste(c[4:length(c)])
					return(strsplit(sentence, " "))
				})
		words <- unlist(words)
		
		# get the unique words and add to result
		terms <- unique(c(terms, unique(words)))
	}
	terms <- sort(terms)
	terms <- terms[!is.na(terms)]
	if(record.secondary.data)
	{	cat("Recording the corpus terms\n")
		dir.create(out.folder,recursive=TRUE,showWarnings=FALSE)
		write.table(x=terms,file=terms.file,row.names=FALSE,col.names=FALSE)
	}
}
cat("\n")

# process each text file
cat("Start extracting the networks\n")
for(i in 1:length(text.files))
#foreach(i=1:length(text.files), .packages="igraph") %dopar%
{	# read the file line-by-line
	text.file <- text.files[i]
	cat("Reading sentences for file ",i,"/",length(text.files)," ",text.file,"\n",sep="")
	sentences <- readLines(paste(in.folder,text.file,sep=""))
	
	# set up output folder
	subfolder <- paste(out.folder,text.file,"/",sep="")
	dir.create(subfolder,recursive=TRUE,showWarnings=FALSE)
	
	# split each line in columns using tab as a separator
	cat("Spliting sentences\n")
	cols <- strsplit(sentences, "\t")
	
	# select the 4th column (sentence)
	words <- lapply(cols,function(col) 
			{	if(length(col)>3)
					paste(col[4:length(col)])
				else
					NA
			})
	
	# possibly remove empty sentences
	idx <- which(sapply(words,is.na))
	if(length(idx)>0)
		words <- words[-idx]

	# split each sentence using space as a separator
	words <- lapply(words,function(w) 
				strsplit(w, " ")[[1]])
	
	# identify and record local terms
	local.terms <- sort(unique(unlist(words)))
	local.terms.indices <- match(local.terms,terms)
	m <- cbind(local.terms.indices,local.terms)
	colnames(m) <- c("GlobalIndex","Term")
	if(record.secondary.data)
		write.table(x=m,file=paste(subfolder,"localterms.txt",sep=""))
	
	# count occurrences
	cat("Counting word occurrences\n")
	#counts <- table(factor(unlist(words),levels=terms))
	counts <- table(factor(unlist(words),levels=local.terms))
	
	# remove one-word sentences
	idx <- which(sapply(words,length)>1)
	words <- words[idx]
	words <- unlist(words)
	
	# build the matrix of adjacent words
	cat("Counting word co-occurrencess\n")
	pairs1 <- lapply(words,function(v) v[1:(length(v)-1)])
	pairs2 <- lapply(words,function(v) v[2:(length(v))])
	pairs <- cbind(unlist(pairs1),unlist(pairs2))
	
	# count co-occurrences
	#co.counts <- table(factor(pairs[,1],levels=terms),factor(pairs[,2],levels=terms))
	#co.counts <- table(factor(pairs[,1],levels=local.terms),factor(pairs[,2],levels=local.terms))
	co.counts <- process.adjacency(mat=pairs, sym=!directed, levels=local.terms)
	
	# record co-occurrence matrix
#	if(record.secondary.data)
	{	cat("Recording the co-occurrence matrix\n")
		write.table(x=co.counts,file=paste(subfolder,"coocurrences.txt",sep=""))
	}
	
	# build the networks
	cat("Building the network\n")
	if(directed)
		mode <- "directed"
	else
		mode <- "undirected"
	g <- graph.adjacency(adjmatrix=co.counts,mode=mode,weighted=TRUE
			#,add.rownames="label" # not necessary (redundant with the 'name' attribute)
		)
	V(g)$Frequency <- counts#graph.strength(g)

	# process centralities
	cat("Processing centralities\n")
	V(g)$Degree <- degree(g)
	V(g)$Betweenness <- betweenness(g)
	V(g)$Closeness <- closeness(g)
	V(g)$Spectral <- evcent(g)$vector
	V(g)$Subgraph <- subgraph.centrality(g)
	
	# process other nodal measures
	cat("Processing other measures\n")
	V(g)$Eccentricity <- eccentricity(g)
	V(g)$Transitivity <- transitivity(graph=g, type="localundirected",isolates="zero")
	
	# detect communities and process related measures
	cat("Processing community-related measures\n")
	coms <- infomap.community(graph=g,modularity=FALSE)
	membr <- membership(coms)
	V(g)$Community <- membr
	V(g)$Embeddedness <- process.embeddedness(g)
	V(g)$GaWithinDeg <- process.ga.withindeg(g)
	V(g)$GaPartCoef <- process.ga.partcoef(g)
	
	# add global measures
	g$Size <- vcount(g)
	g$Density <- graph.density(g)
	g$NormMaxDegree <- max(V(g)$Degree)/vcount(g)
	g$AvrgDegree <- mean(V(g)$Degree)
	g$AvrgBetweenness <- mean(V(g)$Betweenness)
	g$AvrgSpectral <- mean(V(g)$Spectral)
	g$AvrgSubgraph <- mean(V(g)$Subgraph)
	g$AvrgEccentricity <- mean(V(g)$Eccentricity)
	g$AvrgTransitivity <- mean(V(g)$Transitivity)
	g$AvrgEmbeddedness <- mean(V(g)$Embeddedness)
	g$AvrgGaWithinDeg <- mean(V(g)$GaWithinDeg)
	g$AvrgGaPartCoef <- mean(V(g)$GaPartCoef)
	
	# record the network (including all available info)
	if(record.secondary.data)
	{	cat("Recording the network\n")
		write.graph(graph=g,file=paste(subfolder,"wordnetwork.graphml",sep=""),format="graphml")
	}
	
	# record all feature vectors in a single file
	cat("Recording the data\n")
	meas.names <- c("Frequency",
			"Degree","Betweenness","Closeness","Spectral","Subgraph",
			"Eccentricity","Transitivity",
			"Community","Embeddedness","GaWithinDeg","GaPartCoef")
	if(output.full.matrix)
		indices <- local.terms.indices
	else
		indices <- 1:vcount(g)
	data <- matrix(0,nrow=length(indices),ncol=length(meas.names))
	colnames(data) <- meas.names
	for(meas.name in meas.names)
		data[indices,meas.name] <- get.vertex.attribute(graph=g,name=meas.name)
	data <- round(data,digits=4) # limit to 4 decimals
	if(output.full.matrix)
		write.table(x=format(data,scientific=FALSE),file=paste(subfolder,"local.features.txt",sep=""),row.names=FALSE,col.names=TRUE, quote=FALSE, sep="\t")
	else
	{	rownames(data) <- local.terms
		write.table(x=format(data,scientific=FALSE),file=paste(subfolder,"local.features.txt",sep=""),row.names=TRUE,col.names=TRUE, quote=FALSE, sep="\t")
	}

	# plot graph
#	plot(g,vertex.size=4,vertex.label="")

	# record global measures in a separate file
	meas.names <- c("Size","Density",
			"NormMaxDegree","AvrgDegree","AvrgBetweenness","AvrgSpectral","AvrgSubgraph",
			"AvrgEccentricity","AvrgTransitivity",
			"AvrgEmbeddedness","AvrgGaWithinDeg","AvrgGaPartCoef")
	data <- sapply(meas.names,function(meas.name)
				get.graph.attribute(graph=g,name=meas.name))
	names(data) <- meas.names
	data <- round(data,digits=4) # limit to 4 decimals
	write.table(x=format(data,scientific=FALSE),file=paste(subfolder,"global.features.txt",sep=""),row.names=TRUE,col.names=FALSE, quote=FALSE, sep="\t")
}

# stop parallel mode
stopCluster(par.clust)


# process inter-user distances based on the graphs
cat("Start processing the distances\n")
# init result vectors
result <- c()
names1 <- c()
names2 <- c()
# treat each pair of users (i.e. graphs)
for(i in 1:(length(text.files)-1))
{	# set up folder
	text.file1 <- text.files[i]
	name1 <- strsplit(text.file1,"_texts.tsv.out")[[1]]
	subfolder1 <- paste(out.folder,text.file1,"/",sep="")
	cat("..Processing file ",text.file1," (",i,"/",(length(text.files)-1),")\n",sep="")
	
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
	}
}	
# record vectors
df <- data.frame("User1"=names1, "User2"=names2, "Distance"=result)
write.table(x=df,file=paste(out.folder,"distances.txt",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE)

# disable logging
sink()
