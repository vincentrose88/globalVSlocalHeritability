#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=T)

genoSum <- function(x,beta,CAF){
    if(all(x==0,na.rm=T) | all(x==2,na.rm=T)){ #Checking if SNP is monozygotic
        genoWeighted <- NA
    }else{   
        w <- (x-2*CAF)/sqrt(2*CAF*(1-CAF)) 
        genoWeighted <- w*beta
    }
    return(genoWeighted)
}

geno <- read.table(args[1],h=T,as.is=T)

#Remove _Allele from colnames:
colnames(geno)[grep('_',colnames(geno),fixed=T)] <- sapply(colnames(geno)[grep('_',colnames(geno),fixed=T)],function(x){substr(x,1,nchar(x)-2)})
#Read in betas sampled from N(0,0.2)
beta <- read.table(args[2],h=F,as.is=T)

#Read in Freq spectrum
frq <- read.table(args[3],h=T,as.is=T)

betaFrq <- merge(frq[,c('SNP','MAF')],beta,by.x='SNP',by.y='V1')
colnames(betaFrq)[3] <- 'beta'

#Change ':' in chr:pos SNP names to '.' to avoid mismerge with geno
betaFrq[grep(':',betaFrq$SNP,fixed=T),'SNP'] <- sapply(betaFrq[grep(':',betaFrq$SNP,fixed=T),'SNP'],function(x){sub(':','.',x,fixed=T)})

geno.with.beta <- geno #will add stuff to this one
betaFrq.in.geno <- betaFrq[betaFrq$SNP %in% colnames(geno),] #Safety - all should be there

for(i in 1:nrow(betaFrq.in.geno)){
    geno.with.beta[,betaFrq.in.geno[i,'SNP']] <- genoSum(x=geno[,betaFrq.in.geno[i,'SNP']],
                                                         beta=betaFrq[i,'beta'],
                                                         CAF=betaFrq[i,'MAF'])
}

#from 7 and on because the first 6 columns are FID IID PID MID SEX PHE
geno.with.beta$genoBetaSum <- apply(geno.with.beta[,7:(ncol(geno.with.beta)-1)],1,sum,na.rm=T)
t <- geno.with.beta[,c('FID','IID','genoBetaSum')]

write.table(t,args[4],row.names=F,quote=F)


