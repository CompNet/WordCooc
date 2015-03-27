#######################################################
# Functions used to process community-related measures.
# 
# Vincent Labatut 3/3/15
#######################################################
library("igraph")

#######################################################
# Processes the communities of each node neighbors
#######################################################
process.neigh.coms <- function(graph)
{	# process neighborhoods
	neigh <- neighborhood(graph, order=1)
	# process the associated communities
	coms <- lapply(1:vcount(graph),function(u)
			{	nall <- neigh[[u]][neigh[[u]]!=u]
				V(graph)$Community[nall]
			})
	return(coms)	
}

#######################################################
# Processes the embeddedness measure: internal
# degree divided by total degree. The internal
# degree is the number of neighbors in the
# same community.
#
# TODO: could be generalized to handle weights
#######################################################
process.embeddedness <- function(graph)
{	coms <- process.neigh.coms(graph)
	internal.degree <- sapply(1:vcount(graph), function(u)
			{	neigh.coms <- coms[[u]]
				own.com <- V(graph)$Community[u]
				length(which(neigh.coms==own.com))
			})
	result <- internal.degree / V(graph)$Degree
	
	# possibly replace NaN by zeros
	idx <- which(is.nan(result))
	if(length(idx)>0)
		result[idx] <- 0
	
	return(result)
}

###############################################################################
# Processes the community z-score, i.e. the z-score of a given value for each
# node, but processed relatively to the node community (by opposition to the
# whole network.
#
# values: the values to consider.
# membership: communities of the nodes (integer vector).
###############################################################################
process.community.zscores <- function(values, membership)
{	result <- rep(NA,length(values))
	coms <- sort(unique(membership))
	
	for(i in 1:length(coms))
	{	com <- coms[i]
		idx <- which(membership==com)
		result[idx] <- scale(values[idx])
	}
	
	# possibly replace NaN by zeros
	idx <- which(is.nan(result))
	if(length(idx)>0)
		result[idx] <- 0
	
	return(result)
}

#######################################################
# Processes Guimerà & Amaral's Within-degree measure,
# a role measure corresponding to the z-score of the 
# intenal degree.
#
# TODO: could also be generalized to handle weights
#
# graph: the graph to process (community detection is
# 	     supposed to have been processed already).
#######################################################
process.ga.withindeg <- function(graph)
{	coms <- process.neigh.coms(graph)
	internal.degree <- sapply(1:vcount(graph), function(u)
			{	neigh.coms <- coms[[u]]
				own.com <- V(graph)$Community[u]
				length(which(neigh.coms==own.com))
			})
	result <- process.community.zscores(values=internal.degree,V(graph)$Community)
	return(result)
}

#######################################################
# Processes Guimerà & Amaral's Participation coefficient,
# a role measure reflecting the external connectivity,
# i.e. how the node is connected to other communities than 
# its own.
#
# TODO: could also be generalized to handle weights
#
# graph: the graph to process (community detection is
# 	     supposed to have been processed already).
#######################################################
process.ga.partcoef <- function(graph)
{	coms <- process.neigh.coms(graph)
	result <- sapply(1:vcount(graph), function(u)
			{	neigh.coms <- coms[[u]]
				numerator <- (table(neigh.coms))^2
				denominator <- rep((length(neigh.coms))^2,length(numerator))
				res <- 1 - sum(numerator/denominator)
			})
	
	# possibly replace NaN by zeros
	idx <- which(is.nan(result))
	if(length(idx)>0)
		result[idx] <- 0
	
	return(result)
}
