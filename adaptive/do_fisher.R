args<-commandArgs(TRUE)

a1 = args[1]
a2 = args[2]
b1 = args[3]
b2 = args[4]

x = c(a1, a2)
y = c(b1, b2)

z = fisher.test(cbind(x,y))

results = paste(z$p.value, z$estimate, z$conf.int[1:2], sep="\t")

print(results)
