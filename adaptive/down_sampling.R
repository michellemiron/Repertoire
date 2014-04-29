# downsample adaptive data

downsample <- function(tcr, target = 2000000, tnoise = 0.125,  nsim = 10) {
	# target is the number of targeted reads per sample
	# noise is \sigma / \mu for read number
        # nsim: number of simulations

	nc = ncol(tcr)
	nsample = nc - 2
	nclone = nrow(tcr)
	
	results = rep(0, nsim + 1)
	results[1] = jensen_shannon(tcr[1:1000,3], tcr[1:1000,4])

	nreads = rnorm(nsim * 2 , target, tnoise * target)
	p = matrix(ncol=2, nrow=nclone)

	for (i in 1:nsim) {
		for (j in 1:2) {
		    
		    lambda = nreads[(i-1)*2+j] * tcr[,j+2] / sum(tcr[,j+2])
	    	    p[,j] = rpois(nclone, lambda) 
		}

		q = p[p[,1] > 1 | p[,2]>1, ] 
		results[1+i] = jensen_shannon(q[1:1000,1], q[1:1000,2])		
		
	}
#	tcr[,2] = apply(tcr[,3:4], 1, sum)

	cat(results, "\n", sep="\t")

#	return(tcr[tcr$total >=2, ])

}



shannon.entropy <- function(p)
{
	if (min(p) < 0 || sum(p) <= 0)
	   return(NA)
	   p.norm <- p[p>0]/sum(p)
	   -sum(log2(p.norm)*p.norm)
}


jensen_shannon <- function(p, q){
	       ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
 	       # H(X) = \sum x_i * log2(x_i)
#	       p = p[p >0 & q >0]
#	       q = q[p>0 & q>0]
	       p = p / sum(p)
	       q = q / sum(q)
	       Hj = shannon.entropy(0.5 *(p+q)) 
	       Hp = shannon.entropy(p) 
	       Hq = shannon.entropy(q)
	       
	       jsd = Hj - 0.5*(Hp+Hq)
#	       cat(Hj, Hp, Hq, jsd, "\n")
	       return(jsd)
}


# args<-commandArgs(TRUE)
library(getopt)

spec <- matrix(c(
        'input'     , 'i', 1, "character", "input csv or tsv file (required)",
        'type'     , 't', 2, "character", "csv or tsv (optional); default: tsv",
        "noise" , "s", 2, "double", "[0.125] ratio of stddev to mean of reads number",
	"nreads", "n", 2, "double", "[2000000] mean number of reads per sample in each simulation ",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$input)) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}


file = opt$input

itype = "t"
if (!is.null(opt$type)) {
   itype = opt$type
   
}

noise = 0.125
nreads = 2000000
if (!is.null(opt$noise)) {
   noise = opt$noise
}

if(!is.null(opt$nreads)) {
   nreads = opt$nreads			
}


if( itype == "c")  {
    tcr = read.table(file, header=T, sep=",")
} else if (itype == "t") {
  tcr = read.table(file, header=T, sep="\t")			
} else {
  tcr = read.table(file, header=T)			
}

downsample(tcr, nreads , noise)


## write.table(tcr, "", quote=F, sep="\t", row.names = F)




