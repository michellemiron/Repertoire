cloneCal <- function(array, topN = 2000) {
	x = array[array >0 ] / sum(array)
	x = sort(array, decreasing=T)
	
	l = min(length(x), topN )
	x = x[1:l]

	x = x/ sum(x)	
#	print(x[l])
	entropy = sum(x * -1 * log2(x))
	maxentropy = -log2(1/l)
	return(signif(1 - entropy / maxentropy, 2))
}

