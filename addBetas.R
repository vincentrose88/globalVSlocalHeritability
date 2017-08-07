#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=T)

bim <- read.table(args[1],h=F,as.is=T)
snps <- bim$V2
beta <- rnorm(length(snps),0,0.2)
t <- data.frame(snp=snps,beta=beta)
tails <- quantile(beta,c(0.025,0.975))
t <- t[t$beta >= tails[2] | t$beta <= tails[1],]
write.table(t[,1],args[2],row.names=F,col.names=F,quote=F)
write.table(t,paste(args[2],'withBetas',sep='.'),row.names=F,col.names=F,quote=F)
