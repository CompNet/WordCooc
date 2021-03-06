#######################################################
# Processing of the transcription files:
# - Retrieve the text files.
# - Extract the co-occurrence networks.
# - Process a bunch of nodal measures.
# - Record the resulting networks.
# 
# In a file, each line corresponds to a full transcription.
# The script extracts one network for each transcription (line)
# and produce a first file containing the co-occurrence matrix and
# a second one containing the network and the nodal topological measures.
#
# Vincent Labatut 3/2015
#
# setwd("C:/Eclipse/workspaces/Extraction/")
# setwd("~/Eclipse/workspaces/Extraction/")
# source("WordCooc/src/tweets/word-net.R")
#######################################################
library("igraph")

library("compiler")		# enables compilation
enableJIT(3)

source("WordCooc/src/common/com-measures.R")
source("WordCooc/src/common/misc.R")


# Whether or not to record secondary data such as co-occurrence
# networks, list of terms, co-occurrence matrix, etc.
record.secondary.data <- FALSE

# set up in/out folders
in.folder <- "WordCooc/in/clean2/"
out.folder <- "WordCooc/out/nodalmeas/"

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

##################################################################
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
				#table(factor(m[,1],levels=lev),factor(m[,2],levels=lev))
				process.adjacency(mat=m, sym=TRUE, levels=lev)
			})
	
	# record co-occurrence matrices
	if(record.secondary.data)
	{	cat("recording co-occurrence matrices\n")
		sapply(1:length(co.counts), function(i) 
				{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
					dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
					write.table(x=co.counts[[i]],file=paste(sentence.folder,"cooccurrences.txt",sep=""))
				})
	}
	
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
			{	V(g)$Degree <- degree(g)
				V(g)$Betweenness <- betweenness(g)
				V(g)$Closeness <- closeness(g)
				V(g)$Spectral <- evcent(g)$vector
				V(g)$Subgraph <- subgraph.centrality(g)
				return(g)
			})
	
	# process other nodal measures
	cat("Processing other measures\n")
	nets <- lapply(nets,function(g) 
			{	V(g)$Eccentricity <- eccentricity(g)
				V(g)$Transitivity <- transitivity(graph=g, type="localundirected",isolates="zero")
				return(g)
			})
	
	# detect communities and process related measures
	cat("Processing community-related measures\n")
	nets <- lapply(nets,function(g) 
			{	coms <- infomap.community(graph=g,modularity=FALSE)
				membr <- membership(coms)
				V(g)$Community <- membr
				V(g)$Embeddedness <- process.embeddedness(g)
				V(g)$GaWithinDeg <- process.ga.withindeg(g)
				V(g)$GaPartCoef <- process.ga.partcoef(g)
				return(g)
			})
	
	# record the networks (including all available info)
	cat("Recording networks\n")
	nets <- lapply(1:length(nets),function(i) 
			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
				write.graph(graph=nets[[i]],file=paste(sentence.folder,"word-network.graphml",sep=""),format="graphml")
			})
}
