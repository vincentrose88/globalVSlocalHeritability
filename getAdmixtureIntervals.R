#!/usr/local/bin/Rscript

admix <- read.table('admixture.k2.table',h=T,as.is=T)
fam <- read.table('geno/GreenlandSiblings.autosome.geno.0.01.fam',h=F,as.is=T)
fam <- fam[,1:2]

t <- merge(fam,admix,by.x='V1',by.y='id')

for(i in seq(0.5,0.8,0.1)){
    j <- i+0.1
    pop <- t[which(i < t$inuit & t$inuit <= j),1:2]
    write.table(pop,
                paste(i,'to',j,'inuit',sep='.'),
                row.names=F,col.names=F,quote=F)
}

i <- 0.99
j <- 1
pop <- t[which(i < t$inuit & t$inuit <= j),1:2]
write.table(pop,
            paste(i,'to',j,'inuit',sep='.'),
            row.names=F,col.names=F,quote=F)

i <- 0.90
j <- 0.99
pop <- t[which(i < t$inuit & t$inuit <= j),1:2]
write.table(pop,
            paste(i,'to',j,'inuit',sep='.'),
            row.names=F,col.names=F,quote=F)


for(i in seq(0.5,0.7,0.1)){
    j <- i+0.1
    pop <- t[which(i < t$eu & t$eu <= j),1:2]
    write.table(pop,
                paste(i,'to',j,'danes',sep='.'),
                row.names=F,col.names=F,quote=F)
}

#Thresholds
for(i in c(0.99,seq(0.5,1,0.1))){
    pop <- t[which(i < t$inuit),1:2]
    write.table(pop,
                paste(i,'inuit',sep='.'),
                row.names=F,col.names=F,quote=F)
}
#Thresholds
for(i in seq(0.5,0.7,0.1)){
    pop <- t[which(i < t$eu),1:2]
    write.table(pop,
                paste(i,'danes',sep='.'),
                row.names=F,col.names=F,quote=F)
}
