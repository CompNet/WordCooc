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

# TODO we could also consider directed networks (and all the adapted measures)

# Whether or not to record secondary data such as co-occurrence
# networks, list of terms, co-occurrence matrix, etc.
record.secondary.data <- FALSE
# TRUE to record full vectors (including many zeros), i.e. a value
# for each word detected in the corpus, even if it does not appear
# anywhere in the considered user's tweets. FALSE allows focusing
# only on the words actually used.
output.full.matrix <- FALSE

# set up in/out folders
log.file <- "WordCooc/out/log.txt"
#in.folder <- "WordCooc/in/tweets/"
#out.folder <- "WordCooc/out/tweets/"
#in.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/training/out_vincent/"
#out.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/training/out_vincent_features/"
#in.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/test/out_vincent/"
#out.folder <- "/home/imagiweb/works/Replab2014/replab2014_corpus_training_testunlabeled/author_profiling/test/out_vincent_features/"
in.folder <- "~/work/data/training/"
out.folder <- "~/work/data/training_features/"
#in.folder <- "~/work/data/test/"
#out.folder <- "~/work/data/test_features/"

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
	if(record.secondary.data)
	{	cat("Recording the corpus terms\n")
		dir.create(out.folder,recursive=TRUE,showWarnings=FALSE)
		write.table(x=terms,file=terms.file,row.names=FALSE,col.names=FALSE)
	}
}
cat("\n")

# process each text file
cat("Start extracting the networks\n")
#for(i in 1:length(text.files))
foreach(i=1:length(text.files), .packages="igraph") %dopar%
{	# read the file line-by-line
	text.file <- text.files[i]
	cat("Reading sentences for file ",i,"/",length(text.files)," ",text.file,"\n",sep="")
	sentences <- readLines(paste(in.folder,text.file,sep=""))
	
	# set up output folder
	subfolder <- paste(out.folder,text.file,"/",sep="")
	dir.create(subfolder,recursive=TRUE,showWarnings=FALSE)
	
	# remove empty sentences
	sentences <- sentences[nchar(sentences,type="bytes")>0]
	
	# split each line in columns using tab as a separator
	cat("Spliting sentences\n")
	cols <- strsplit(sentences, "\t")
	
	# split each sentence using space as a separator
	words <- lapply(cols,function(c) 
			{	sentence <- paste(c[4:length(c)])
				return(strsplit(sentence, " ")[[1]])
			})
	
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
	co.counts <- process.adjacency(mat=pairs, sym=TRUE, levels=local.terms)
	
	# record co-occurrence matrix
#	if(record.secondary.data)
	{	cat("Recording the co-occurrence matrix\n")
		write.table(x=co.counts,file=paste(subfolder,"coocurrences.txt",sep=""))
	}
	
	# build the networks
	cat("Building the network\n")
	g <- graph.adjacency(adjmatrix=co.counts,mode="undirected",weighted=TRUE
			#,add.rownames="label" # not necessary (redundant with the 'name' attribute)
		)
	V(g)$frequency <- counts#graph.strength(g)

	# process centralities
	cat("Processing centralities\n")
	V(g)$degree <- degree(g)
	V(g)$betweenness <- betweenness(g)
	V(g)$closeness <- closeness(g)
	V(g)$spectral <- evcent(g)$vector
	V(g)$subgraph <- subgraph.centrality(g)
	
	# process other nodal measures
	cat("Processing other measures\n")
	V(g)$eccentricity <- eccentricity(g)
	V(g)$transitivity <- transitivity(graph=g, type="localundirected",isolates="zero")
	
	# detect communities and process related measures
	cat("Processing community-related measures\n")
	coms <- infomap.community(graph=g,modularity=FALSE)
	membr <- membership(coms)
	V(g)$community <- membr
	V(g)$embeddeness <- process.embeddedness(g)
	V(g)$gawithindeg <- process.ga.withindeg(g)
	V(g)$gapartcoef <- process.ga.partcoef(g)
	
	# record the network (including all available info)
	if(record.secondary.data)
	{	cat("Recording the network\n")
		write.graph(graph=g,file=paste(subfolder,"wordnetwork.graphml",sep=""),format="graphml")
	}
	
	# record all feature vectors in a single file
	cat("Recording the data\n")
	meas.names <- c("frequency",
			"degree","betweenness","closeness","spectral","subgraph",
			"eccentricity","transitivity",
			"community","embeddeness","gawithindeg","gapartcoef")
	if(output.full.matrix)
		indices <- local.terms.indices
	else
		indices <- 1:vcount(g)
	data <- matrix(0,nrow=length(indices),ncol=length(meas.names))
	colnames(data) <- meas.names
	for(meas.name in meas.names)
		data[indices,meas.name] <- get.vertex.attribute(graph=g,name=meas.name)
	data <- round(data,digits=4)
	if(output.full.matrix)
		write.table(x=format(data,scientific=FALSE),file=paste(subfolder,"features.txt",sep=""),row.names=FALSE,col.names=TRUE, quote=FALSE, sep="\t")
	else
	{	rownames(data) <- local.terms
		write.table(x=format(data,scientific=FALSE),file=paste(subfolder,"features.txt",sep=""),row.names=TRUE,col.names=TRUE, quote=FALSE, sep="\t")
	}

	# plot graph
#	plot(g,vertex.size=4,vertex.label="")
}
stopCluster(par.clust)

# disable logging
sink()
