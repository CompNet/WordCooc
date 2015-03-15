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

source("WordCooc/com-measures.R")


# set up in/out folders
#in.folder <- "WordCooc/in/tweets/"
in.folder <- "D:/Users/Vincent/Documents/Travail/training-out_vincent/"
#in.folder <- "D:/Users/Vincent/Documents/Travail/test-out_vincent/"
out.folder <- "WordCooc/out/tweets/training/"
#out.folder <- "WordCooc/out/tweets/test/"

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

# build the list of terms for the whole collection
terms.file <- paste(out.folder,"terms.txt",sep="")
if(file.exists(terms.file))
{	cat("Loading the corpus terms\n")
	terms <- as.matrix(read.table(terms.file))	
}else
{	cat("Identifying the corpus terms\n")
	terms <- c()
	for(text.file in text.files)
	{	# read the file
		cat("Reading sentences for",text.file,"\n")
		sentences <- readLines(paste(in.folder,text.file,sep=""))
		
		# remove empty sentences
		sentences <- sentences[nchar(sentences)>0]
		
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
	cat("Recording the corpus terms\n")
	dir.create(out.folder,recursive=TRUE,showWarnings=FALSE)
	write.table(x=terms,file=terms.file,row.names=FALSE,col.names=FALSE)
}
cat("\n")

# process each text file
cat("Start extracting the networks\n")
for(text.file in text.files)
{	# read the file line-by-line
	cat("Reading sentences for",text.file,"\n")
	sentences <- readLines(paste(in.folder,text.file,sep=""))
	
	# set up output folder
	subfolder <- paste(out.folder,text.file,"/",sep="")
	dir.create(subfolder,recursive=TRUE,showWarnings=FALSE)
	
	# remove empty sentences
	sentences <- sentences[nchar(sentences)>0]
	
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
	write.table(x=m,file=paste(subfolder,"localterms.txt",sep=""))
	
	# count occurrences
	cat("Counting word occurrences\n")
#	counts <- table(factor(unlist(words),levels=terms))
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
#	co.counts <- table(factor(pairs[,1],levels=terms),factor(pairs[,2],levels=terms))
	co.counts <- table(factor(pairs[,1],levels=local.terms),factor(pairs[,2],levels=local.terms))
	
	# record co-occurrence matrix
	cat("Recording the co-occurrence matrix\n")
	write.table(x=co.counts,file=paste(subfolder,"coocurrences.txt",sep=""))
	
	# build the networks
	cat("Building the network\n")
	g <- graph.adjacency(adjmatrix=co.counts,mode="undirected",weighted=TRUE
#			,add.rownames="label" # not necessary (redundant with the 'name' attribute)
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
	cat("Recording the network\n")
	write.graph(graph=g,file=paste(subfolder,"wordnetwork.graphml",sep=""),format="graphml")
	
	# record the full vectors in a single file
	cat("Recording the data\n")
	meas.names <- c("frequency",
			"degree","betweenness","closeness","spectral","subgraph",
			"eccentricity","transitivity",
			"community","embeddeness","gawithindeg","gapartcoef")
	data <- matrix(0,ncol=length(meas.names),nrow=length(terms))
	colnames(data) <- meas.names
	for(meas.name in meas.names)
		data[local.terms.indices,meas.name] <- get.vertex.attribute(graph=g,name=meas.name)
	write.table(x=data,file=paste(subfolder,"features.txt",sep=""),row.names=FALSE,col.names=TRUE)
	
#	con <- file(paste(subfolder,"features.txt",sep=""),open="wt")
#	record.vector <- function(vect, con, sz, idx)
#	{	v <- rep(0,sz)
#		v[idx] <- vect
#		writeLines(text=v,con)
#	}
#	record.vector(vect=V(g)$frequency, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$degree, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$betweenness, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$closeness, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$spectral, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$subgraph, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$eccentricity, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$transitivity, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$community, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$embeddeness, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$gawithindeg, con=con, sz=length(terms), idx=local.terms.indices)
#	record.vector(vect=V(g)$gapartcoef, con=con, sz=length(terms), idx=local.terms.indices)
#	close(con)

	# plot graph
#	plot(g,vertex.size=4,vertex.label="")
}
