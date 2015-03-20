#######################################################
# Processing of the transcription files:
# - Retrieve the text files.
# - Retrieve the topic map.
# - Extract the co-occurrence networks.
# - Replace words by the corresponding topic.
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
# source("WordCooc/extract-net-transcriptions-topics.R")
#######################################################
library("igraph")

source("WordCooc/misc.R")


# States whether word order should be considered or not.
# If TRUE, then the sequence word1 word2 leads to a link
# from word1 to word2, but not in the other direction. 
# If FALSE, both sequences word1 word2 and word2 word1 are
# considered similarly.
directed <- TRUE 

# States how to consider separators. Separators are words belonging 
# to no topic at all. If "ignore", they are just removed from the
# text before process, and do not appear in the final network.
# If "explicit", they appear as an additional topic node in the
# network. If "implicit", they are considered as separators (two
# topic words separated by a no-topic word are not considered as
# adjacent), but do not appear in the topic network.
separator <- "ignore" # ignore explicit implicit 

# set up in/out folders
in.folder <- "WordCooc/in/clean2/"
in.map <- "WordCooc/in/discriminativeWordsListForEachTheme_200.txt"
out.folder <- "WordCooc/out/topics/"

# define file prefix (for the generated files)
prefix <- paste(
		"separator=",separator,".",
		"directed=",tolower(directed),".",
		sep="")

# get the topic map
topic.map <- as.matrix(read.table(in.map))
topic.names <- sort(unique(topic.map[,2]))
topic.map <- rbind(topic.map,c("NA","SEP")) # add the separator symbol
if(separator=="explicit")
	topic.names <- sort(unique(topic.map[,2]))

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

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
	
	# replace each word by its topic
	cat("Replacing words by corresponding topics\n")
	topics <- lapply(words,function(v) 
			{	idx <- match(v,topic.map[,1])	# get the row of the word in the map
				idx[is.na(idx)] <- nrow(topic.map) # for words not contained in the map, use the last row (=SEP)
				return(topic.map[idx,2])
			})
	
	# possibly remove all no-topic words
	if(separator=="ignore")
	{	topics <- lapply(topics,function(v) 
				{	idx <- which(v=="SEP")
					v <- v[-idx]
				})
		# possibly update indices for empty sentences
		idx.rmd <- which(sapply(topics,length)==0)
		if(length(idx.rmd)>0)
		{	idx.kpt <- idx.kpt[-idx.rmd]
			topics <- topics[-idx.rmd]
		}
	}
	
	# build the matrices of adjacent words
	cat("Counting topic co-occurrencess\n")
	pairs <- lapply(topics,function(v) 
				cbind(v[1:(length(v)-1)],v[2:(length(v))]))
	
	# build the adjacency matrices
	co.counts <- lapply(pairs, function(m) 
				process.adjacency(mat=m, sym=!directed, levels=topic.names))
		
	# record co-occurrence matrices (only the non-redundant part)
	cat("recording co-occurrence matrices as vectors\n")
	sapply(1:length(co.counts), function(i) 
			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
				m <- co.counts[[i]]
				if(!directed)							# if the graph is undirected, the matrix is symmetrical
					m[lower.tri(m,diag=FALSE)] <- NA	# put NA in the lower triangle of the matrix
				m <- as.data.frame(as.table(m))  		# turn into a 3-column table
				m <- na.omit(m)							# possibly remove NAs
				data <- m[,3]
				names(data) <- paste(m[,1],m[,2],sep="-")
				write.table(x=data,file=paste(sentence.folder,prefix,"coocurrences.txt",sep=""),col.names=FALSE,quote=FALSE)
			})

	# build the networks
#	cat("Building networks\n")
#	nets <- lapply(1:length(co.counts),function(i) 
#				g <- graph.adjacency(adjmatrix=co.counts[[i]],mode="undirected",weighted=TRUE))

	# record the networks (including all available info)
#	cat("Recording networks\n")
#	nets <- lapply(1:length(nets),function(i) 
#			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
#				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
#				write.graph(graph=nets[[i]],file=paste(sentence.folder,prefix,"topic-network.graphml",sep=""),format="graphml")
#			})
}
