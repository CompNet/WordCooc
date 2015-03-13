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
# setwd("C:/Eclipse/workspaces/Extraction/")
# setwd("~/Eclipse/workspaces/Extraction/")
# source("WordCooc/extract-net-tweets.R")
#######################################################
library("igraph")

source("WordCooc/com-measures.R")


# set up in/out folders
in.folder <- "WordCooc/in/tweets/"
out.folder <- "WordCooc/out/tweets/"

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

# build the list of terms for the whole collection
terms.file <- paste(out.folder,"terms.txt",sep="")
if(file.exits(terms.file))
{	cat("Loading the corpus terms\n")
	terms <- as.matrix(read.table(terms.file))	
}else
{	cat("Identifying the corpus terms\n")
	terms <- c()
	for(text.file in text.files)
	{	# read the file
		cat("Reading sentences\n")
		sentences <- readLines(paste(in.folder,text.file,sep=""))
		
		# remove empty sentences
		sentences <- sentences[nchar(sentences)>0]
		
		# split each line using space as a separator
		cat("Spliting sentences\n")
		words <- strsplit(sentences, " ")
		words <- unlist(words)
		
		# get the unique words and add to result
		terms <- unique(c(terms, unique(words)))
	}
	terms <- sort(terms)
	cat("Recording the corpus terms\n")
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
	
	# remove empty sentences
	sentences <- sentences[nchar(sentences)>0]
	
	# split each line using space as a separator
	cat("Spliting sentences\n")
	words <- strsplit(sentences, " ")
	
	# count occurrences
	cat("Counting word occurrences\n")
	counts <- table(factor(unlist(words),levels=terms))
	
	# build the matrix of adjacent words
	cat("Counting word co-occurrencess\n")
	pairs1 <- lapply(words,function(v) v[1:(length(v)-1)])
	pairs2 <- lapply(words,function(v) v[2:(length(v))])
	pairs <- cbind(unlist(pairs1),unlist(pairs2))
	
	# count co-occurrences
	co.counts <- table(factor(pairs[,1],levels=terms),factor(pairs[,2],levels=terms))
	
	# record co-occurrence matrix
	cat("Recording the co-occurrence matrix\n")
	dir.create(subfolder,recursive=TRUE,showWarnings=FALSE)
	write.table(x=co.counts,file=paste(subfolder,"coocurrences.txt",sep=""))
	
	# build the networks
	cat("Building the network\n")
	g <- graph.adjacency(adjmatrix=co.counts,mode="undirected",weighted=TRUE
#			,add.rownames="label" # not necessary (redundant with the 'name' attribute)
		)
	V(g)$frequency <- counts[[i]]#graph.strength(g)

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
	dir.create(subfolder,recursive=TRUE,showWarnings=FALSE)
	write.graph(graph=g,file=paste(subfolder,"wordnetwork.graphml",sep=""),format="graphml")
}
