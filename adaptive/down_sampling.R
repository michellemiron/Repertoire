# downsample adaptive data

downsample <- function(tcr, target = 2000000) {
	# target is the number of targeted reads per sample
	total = sum(tcr$total)
	nc = ncol(tcr)
	nsample = nc - 2
	nclone = nrow(tcr)
	
	ratio = target / total * nsample
	
	# Bernoulli 
	for (i in 1:nclone) {
	       	ratioA = ratio * tcr[i,3:nc] / sum(tcr[i,3:nc])
	       	for (j in 1:nsample) {
		    tcr[i, j + 2] = round(rbinom(1, tcr$total[i], ratioA[j]))
		}
		tcr[i, 2] = sum(tcr[i, 3:nc])
		
	}
	
	# normalize
	for (j in 3:nc) {
		tcr[,j] = tcr[,j] / sum(tcr[,j]) * 100
	}
	return(tcr)

}





