clonality <- function(array) {
	x = normalize(array)
	l = length(x)
	entropy = sum(x * -1 * log2(x))
	maxentropy = -log2(1/l)
	return(signif(1 - entropy / maxentropy, 2))
}


simpsonIndex <- function(data) {
	a = normalize(data) # data[data >0]
	si = sum(a**2)
	return(si)

}

normalize <- function(data) {
	data = data[data > 0]
	data = data / sum(data)
	
	return(data)
}
