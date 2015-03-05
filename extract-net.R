#######################################################
# Main processing:
# - Retrieve the text files.
# - Extract the co-occurrence networks.
# - Process a bunch of nodal measures.
# - Record the resulting networks.
# 
# Vincent Labatut 3/2015
#######################################################
library("igraph")

source("WordCooc/com-measures.R")


# set up in/out folders
in.folder <- "WordCooc/in/clean/"
out.folder <- "WordCooc/out/"

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)

# process each text file
for(text.file in text.files)
{	# read the file line-by-line
	cat("Reading sentences\n")
	sentences <- readLines(paste(in.folder,text.file,sep=""))
	
	# set up output folder
	subfolder <- paste(out.folder,text.file,"/",sep="")
	
	# remove empty sentences
	idx.kpt <- 1:length(sentences)
	idx.rmd <- which(nchar(sentences)==0)
	if(length(idx.rmd)>0)
	{	idx.kpt <- idx.kpt[-idx.rmd]
		sentences <- sentences[-idx.rmd]
	}
	
	# split each line using space as a separator
	cat("Spliting sentences\n")
	words <- strsplit(sentences, " ")
	
	# count occurrences
	cat("Counting word occurrences\n")
	counts <- lapply(words,function(v) 
			{	lev = sort(unique(v))
				table(factor(v,levels=lev))
			})
	
	# build the matrices of adjacent words
	cat("Counting word co-occurrencess\n")
	pairs <- lapply(words,function(v) 
				cbind(v[1:(length(v)-1)],v[2:(length(v))]))
	
	# count co-occurrences
	co.counts <- lapply(pairs,function(m) 
			{	lev = sort(unique(c(m)))
				table(factor(m[,1],levels=lev),factor(m[,2],levels=lev))
			})
	
	# record co-occurrence matrices
	cat("Cording co-occurrence matrices\n")
	sapply(1:length(co.counts), function(i) 
			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
				write.table(x=co.counts[[i]],file=paste(sentence.folder,"coocurrences.txt",sep=""))
			})
	
	# build the networks
	cat("Building networks\n")
	nets <- lapply(1:length(co.counts),function(i) 
			{	g <- graph.adjacency(adjmatrix=co.counts[[i]],mode="undirected",weighted=TRUE
#					,add.rownames="label" # not necessary (redundant with the 'name' attribute)
				)
				V(g)$frequency <- counts[[i]]#graph.strength(g)
				return(g)
			})

	# process centralities
	cat("Processing centralities\n")
	nets <- lapply(nets,function(g) 
			{	V(g)$degree <- degree(g)
				V(g)$betweenness <- betweenness(g)
				V(g)$closeness <- closeness(g)
				V(g)$spectral <- evcent(g)$vector
				V(g)$subgraph <- subgraph.centrality(g)
				return(g)
			})
	
	# process other nodal measures
	cat("Processing other measures\n")
	nets <- lapply(nets,function(g) 
			{	V(g)$eccentricity <- eccentricity(g)
				V(g)$transitivity <- transitivity(graph=g, type="localundirected",isolates="zero")
				return(g)
			})
	
	# detect communities and process related measures
	cat("Processing community-related measures\n")
	nets <- lapply(nets,function(g) 
			{	coms <- infomap.community(graph=g,modularity=FALSE)
				membr <- membership(coms)
				V(g)$community <- membr
				V(g)$embeddeness <- process.embeddedness(g)
				V(g)$gawithindeg <- process.ga.withindeg(g)
				V(g)$gapartcoef <- process.ga.partcoef(g)
				return(g)
			})
	
	# record the networks (including all available info)
	cat("Recording networks\n")
	nets <- lapply(1:length(nets),function(i) 
			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
				write.graph(graph=nets[[i]],file=paste(sentence.folder,"wordnetwork.graphml",sep=""),format="graphml")
			})
}
