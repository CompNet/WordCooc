#######################################################
# Various functions used in the other scripts.
# 
# Vincent Labatut 3/3/15
#######################################################


#######################################################
# Processes an adjacency table, possibly making it
# symmetrical.
#
# mat: a 2-column matrix (=links).
# sym: whether the adjacency matrix should be
# 	   symmetrical or not.
# levels: the levels fetched to the table function.
#######################################################
process.adjacency <- function(mat, sym, levels)
{	res <- as.matrix(table(factor(mat[,1],levels=levels),factor(mat[,2],levels=levels)))
	if(sym)
	{	res2 <- t(res)
		diag(res2) <- 0
		res <- res + res2
	}
	return(res)
}
