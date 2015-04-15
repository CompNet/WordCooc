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
# source("WordCooc/src/tweets/topic-net.R")
#######################################################
library("igraph")

library("compiler")		# enables compilation
enableJIT(3)

source("WordCooc/src/common/misc.R")


# States whether word order should be considered or not.
# If TRUE, then the sequence word1 word2 leads to a link
# from word1 to word2, but not in the other direction. 
# If FALSE, both sequences word1 word2 and word2 word1 are
# considered similarly.
directed <- FALSE

# States how to consider separators. Separators are words belonging 
# to no topic at all. If "ignore", they are just removed from the
# text before process, and do not appear in the final network.
# If "explicit", they appear as an additional topic node in the
# network. If "implicit", they are considered as separators (two
# topic words separated by a no-topic word are not considered as
# adjacent), but do not appear in the topic network.
separator <- "ignore" # ignore explicit implicit

# Whether or not collapsing consecutive words of the same topic,
# i.e. merging them and representing them by a single occurrence
# of the said topic.
collapsed <- TRUE

# Threshold used to define the topics: minimal relevance value 
# a word must have to be considered as characterisitc of a topic.
threshold <- 0

# set up in/out folders
in.folder <- "WordCooc/in/clean2/"
#in.folder <- "WordCooc/in/test/"
out.folder <- "WordCooc/out/topics/"

# set up topic map
#in.map <- "WordCooc/in/discriminativeWordsListForEachTheme_200.txt"
#in.map <- "WordCooc/in/gini_200.txt"
in.map <- "WordCooc/in/gini_500.txt"



# define file prefix (for the generated files)
prefix <- paste(
		"separator=",separator,".",
		"directed=",tolower(directed),".",
		"collapsed=",tolower(collapsed),".",
		sep="")

# get the topic map
topic.table <- read.table(in.map)
idx <- 1:nrow(topic.table)
if(ncol(topic.table)>2)
	idx <- which(topic.table[,3]>=threshold)
topic.map <- as.matrix(topic.table[idx,1:2])
topic.names <- sort(unique(topic.map[,2]))
topic.map <- rbind(topic.map,c("NA","SEP")) # add the separator symbol
if(separator=="explicit")
	topic.names <- sort(unique(topic.map[,2]))

# set up levels for triplets
trilev <- expand.grid(topic.names,topic.names,topic.names)
if(directed)
{	trilev <- apply(trilev, 1, function(v) paste(rev(v),collapse="-"))
}else
{	trilev <- apply(trilev, 1, function(v) paste(sort(v),collapse="-"))
	trilev <- sort(unique(trilev))
}		

# get text files
text.files <- list.files(path=in.folder,full.names=FALSE,no..=TRUE)
print(text.files)

##################################################################
# process each text file
for(text.file in text.files)
{	####### input/output
	# read the file line-by-line
	cat("Reading sentences\n")
	sentences <- readLines(paste(in.folder,text.file,sep=""))
	
	# set up output folder
	subfolder <- paste(out.folder,text.file,"/",sep="")
	
	####### pre-processing
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
					if(length(idx)>0)
						v <- v[-idx]
					return(v)
				})
		# possibly update indices for empty sentences
		idx.rmd <- which(sapply(topics,length)==0)
		if(length(idx.rmd)>0)
		{	idx.kpt <- idx.kpt[-idx.rmd]
			topics <- topics[-idx.rmd]
		}
	}
	
	# possibly merge similar consecutive topics
	if(collapsed)
	{	topics <- lapply(topics,function(v)
					rle(v)$values)
	}
	
	####### 1 topic
#	# process the frequency of topics
#	cat("Counting topics\n")
#	freq <- lapply(topics,function(v) 
#				as.matrix(table(factor(v,levels=topic.names))))
#	
#	# record topic frequencies
#	cat("Recording topic frequencies\n")
#	sapply(1:length(freq), function(i) 
#		{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
#			dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
#			write.table(x=freq[[i]],file=paste(sentence.folder,"separator=",separator,".freq.txt",sep=""),col.names=FALSE,quote=FALSE)
#		})
	
	####### 2 topics
#	# build the matrices of adjacent words
#	cat("Counting pair co-occurrencess\n")
#	pairs <- lapply(topics,function(v) 
#				cbind(v[1:(length(v)-1)],v[2:(length(v))]))
#	
#	# build the adjacency matrices
#	co.counts <- lapply(pairs, function(m) 
#				process.adjacency(mat=m, sym=!directed, levels=topic.names))
#	
#	# record co-occurrence matrices (only the non-redundant part)
#	cat("Recording co-occurrence matrices as vectors\n")
#	sapply(1:length(co.counts), function(i) 
#			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
#				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
#				m <- co.counts[[i]]
#				if(!directed)							# if the graph is undirected, the matrix is symmetrical
#					m[lower.tri(m,diag=FALSE)] <- NA	# put NA in the lower triangle of the matrix
#				m <- as.data.frame(as.table(m))  		# turn into a 3-column table
#				m <- na.omit(m)							# possibly remove NAs
#				data <- m[,3]
#				names(data) <- paste(m[,1],m[,2],sep="-")
#				write.table(x=data,file=paste(sentence.folder,prefix,"pairs.txt",sep=""),col.names=FALSE,quote=FALSE)
#			})
#	
#	# build the networks
#	cat("Building networks\n")
#	nets <- lapply(1:length(co.counts),function(i) 
#				g <- graph.adjacency(adjmatrix=co.counts[[i]],mode="undirected",weighted=TRUE))
#
#	# record the networks (including all available info)
#	cat("Recording networks\n")
#	lapply(1:length(nets),function(i) 
#			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
#				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
#				write.graph(graph=nets[[i]],file=paste(sentence.folder,prefix,"network.graphml",sep=""),format="graphml")
#			})

	# networks of higher orders
	cat("Processing higher-order networks\n")
	orders.mat <- lapply(topics, function(top)
	{	temp.list <- lapply(0:3, function(k)
		{	# remove certain topics in the original vectors
			partial.topics <- lapply(0:k, function(i)
				if(length(top)>(k+1))
					top[seq(1+i,length(top),k+1)]
				else
					c())
			
			# process the matrices of pairs of (pseudo)consecutive topics
			partial.pairs <- lapply(partial.topics,function(v)
				if(length(v)>1)
					cbind(v[1:(length(v)-1)],v[2:(length(v))])
				else
					c())
			
			# process the corresponding cooccurrence matrices 
			partial.co.counts <- lapply(partial.pairs, function(m)
				process.adjacency(mat=m, sym=!directed, levels=topic.names))
			
			# sum the list of matrices to get a single total matrix
			m <- Reduce('+',partial.co.counts)
			
			# possibly keep only the upper triangle of the total matrix, and linearize it
			if(!directed)							# if the graph is undirected, the matrix is symmetrical
				m[lower.tri(m,diag=FALSE)] <- NA	# put NA in the lower triangle of the matrix
			m <- as.data.frame(as.table(m))  		# turn into a 3-column table
			m <- na.omit(m)							# possibly remove NAs
			data <- m[,3]
			names(data) <- paste(m[,1],m[,2],sep="-")
			return(data)	
		})
		df <- data.frame(temp.list)
		colnames(df) <- c("Order1","Order2","Order3","Order4")
		return(df)
	})
	# record the resulting linearized matrices
	cat("Recording co-occurrence matrices as vectors\n")
	sapply(1:length(orders.mat), function(i) 
		{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
			dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
			write.table(x=orders.mat[[i]],file=paste(sentence.folder,prefix,"orders.txt",sep=""),col.names=FALSE,quote=FALSE)
		})
	# process correlations between columns and average over all texts
	cat("Process correlations:\n")
	correlation.matrices <- lapply(orders.mat, cor)
	m <- Reduce('+',correlation.matrices) / length(correlation.matrices)
	print(m)
	
	####### 3 topics
#	# matrix of triplets
#	cat("Counting triplet co-occurrencess\n")
#	triplets <- lapply(topics,function(v) 
#				cbind(v[1:(length(v)-2)],v[2:(length(v)-1)],v[3:(length(v))]))
#
#	# counting triplets
#	merged <- lapply(triplets, function(m) 
#				apply(m, 1, function(v) 
#						{	if(!directed)
#								v <- sort(v)
#							paste(v,collapse="-")
#						})
#				)
#	coco.counts <- lapply(merged, function(v) 
#				table(factor(v,levels=trilev)))
#	
#	# recording triplet counts
#	lapply(1:length(coco.counts), function(i)
#			{	sentence.folder <- paste(subfolder,idx.kpt[i],"/",sep="")
#				dir.create(sentence.folder,recursive=TRUE,showWarnings=FALSE)
#				m <- coco.counts[[i]]
#				write.table(x=m,file=paste(sentence.folder,prefix,"triplets.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE)
#			})
}
