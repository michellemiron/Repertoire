## make scatter plots from pair-wise read counts table

makeTable <- function(file) {
	t = read.table(file, header=T, sep="\t")
	subcol = c(1,3)
	t = t[,subcol]
	t.name = sub(".tsv", "", file)
	t.name = sub("CD4", "", t.name)
	t.name = sub("SP", " SP", t.name)
	t.name = sub("ILN", " ILN", t.name)
	t.name = sub("LLN", " LLN", t.name)
	t.name = paste("D", t.name, sep ="")
	colnames(t) = c("nucleotide", t.name)
	return(t)

}

shannon.entropy <- function(p)
{
        if (min(p) < 0 || sum(p) <= 0)
                return(NA)
        p.norm <- p[p>0]/sum(p)
        -sum(log2(p.norm)*p.norm)
}


jensen_shannon <- function(x){
		p = x[,1]
		q = x[,2]
        ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
        # H(X) = \sum x_i * log2(x_i)
#       p = p[p >0 & q >0]
#       q = q[p>0 & q>0]
        p = p / sum(p)
        q = q / sum(q)
        Hj = shannon.entropy(0.5 *(p+q)) 
        Hp = shannon.entropy(p) 
        Hq = shannon.entropy(q)
        
        jsd = Hj - 0.5*(Hp+Hq)
#       cat(Hj, Hp, Hq, jsd, "\n")
        return(jsd)
}


plotOnePair <- function(data, cols ) {
	r = range(c(0.5,max(max(data[,2]), max(data[,3], max(data[,4])))))

	plot(data[,cols]+0.5, log="xy", axes = F,  frame.plot=TRUE, xlim=r, ylim=r, col='#2171B5')
	axis(1, at=c(0.5, 10, 100, 1000, 10000), labels=c("0",  "10", "100", "1000", "10000"))
	axis(2, at=c(0.5, 10, 100, 1000, 10000), labels=c("0",  "10", "100", "1000", "10000"))
	
	axis(3, at=c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)*sum(data[,cols[1]]), labels=c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
	axis(4, at=c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)*sum(data[,cols[2]]), labels=c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
	jsd = jensen_shannon(data[,cols])
	# title(main=paste("JSD = ", round(jsd,4), sep=""), font.main=1)
	legend("top", paste("JSD = ", round(jsd, 2), sep="") , bty="n")
}

plotOneDonor <- function(donor) {
	f1 = paste(donor, "SPCD4.tsv", sep="")
	f2 = paste(donor, "LLNCD4.tsv", sep="")
	f3 = paste(donor, "ILNCD4.tsv", sep="")
	
	s1 = makeTable(f1)
	s2 = makeTable(f2)
	s3 = makeTable(f3)
	s1.col = colnames(s1)
	s2.col = colnames(s2)
	s3.col = colnames(s3)
	
	s1 = cbind(s1[1], s1[2], rep(0, nrow(s1)), rep(0, nrow(s1)))
	s2 = cbind(s2[1], rep(0, nrow(s2)), s2[2], rep(0, nrow(s2)))
	s3 = cbind(s3[1], rep(0, nrow(s3)), rep(0, nrow(s3)), s3[2])
	
	colnames(s1) = c(s1.col[1], s1.col[2], s2.col[2], s3.col[2])
	colnames(s2) = colnames(s1)
	colnames(s3) = colnames(s1)
	
	comb = aggregate(. ~ nucleotide, data = rbind(s1, s2, s3), sum)
	
	
	
	sa = c(2,3)
	plotOnePair(comb, c(2,3))
	plotOnePair(comb, c(2,4))
	plotOnePair(comb, c(3,4))

}

## pdf(paste(paste(colnames(comb)[2:3], collapse = '_'), "pdf", sep="."), width = 15, height= 25 )


pdf("CD4.pdf", width=12, height=21)
par(mfrow=c(5,3))
par(oma=c(1,1,1,1))
par(mar=c(5,4,4,4))
#par(mgp=c(2,1,1))

for (d in c("72", "73", "76", "79", "86")) {
#for (d in c("72")) {
	plotOneDonor(d)
}


dev.off()
