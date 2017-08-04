#!/usr/bin/Rscript
localH <- NULL
for(h in c(0.4,0.6,0.8)){
    for(sim in 1:100){
#       sim='selected' #'sameGeno'#
        phe <- read.table(
            paste('trunc.sim/sim',sim,'h',h,'gctaSim.phen',sep='.'),
            h=F,as.is=T)
        colnames(phe) <- c('FID','IID','phe')
        
        all <- read.table(
            paste('trunc.sim/sim',sim,'all.genoSum',sep='.')
            ,h=T,as.is=T)
        
        localMerge <- merge(all,phe,by=c('FID','IID'))
        tmpLocal <- c(h,var(localMerge$genoBetaSum)/var(localMerge$phe),NA,'all',nrow(localMerge),sim)
        names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')

        localH <- rbind(localH,tmpLocal)

        #inuits
        for(i in seq(0.5,0.8,0.1)){
            j <- i+0.1
            admix <- paste(i,'to',j,sep='.')
            pop <- 'inuit' 
           localGenoSum <- read.table(
                paste('trunc.sim/sim',sim,admix,pop,'genoSum',sep='.')
                ,h=T,as.is=T)
            
            localMerge <- merge(localGenoSum,phe,by=c('FID','IID'))
            tmpLocal <- c(h,var(localGenoSum$genoBetaSum)/var(localMerge$phe),admix,pop,nrow(localMerge),sim)
            names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')
            
            localH <- rbind(localH,tmpLocal)
            
       }
        #inuits special
        for(admix in c('0.9.to.0.99','0.99.to.1')){
            pop <- 'inuit'
            localGenoSum <- read.table(
                paste('trunc.sim/sim',sim,admix,pop,'genoSum',sep='.')
                ,h=T,as.is=T)
            
            localMerge <- merge(localGenoSum,phe,by=c('FID','IID'))
            tmpLocal <- c(h,var(localGenoSum$genoBetaSum)/var(localMerge$phe),admix,pop,nrow(localMerge),sim)
            names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')
            
            localH <- rbind(localH,tmpLocal)
            
       }
        #Danes
        for(i in seq(0.5,0.7,0.1)){
            j <- i+0.1
            admix <- paste(i,'to',j,sep='.')
            pop <- 'danes'
            localGenoSum <- read.table(
                paste('trunc.sim/sim',sim,admix,pop,'genoSum',sep='.')
                ,h=T,as.is=T)
            
            localMerge <- merge(localGenoSum,phe,by=c('FID','IID'))
            tmpLocal <- c(h,var(localGenoSum$genoBetaSum)/var(localMerge$phe),admix,pop,nrow(localMerge),sim)
            names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')
            
            localH <- rbind(localH,tmpLocal)                
        }
        
    }
}


#########
# thresholds

for(h in c(0.4,0.6,0.8)){
    for(sim in 1:100){   
        phe <- read.table(
            paste('trunc.sim/sim',sim,'h',h,'gctaSim.phen',sep='.'),
            h=F,as.is=T)
        colnames(phe) <- c('FID','IID','phe')

        #Special case for super inuit (90% and 99%, for which there is no danes equilivant. 
        for(admix in c('0.8','0.9','0.99')){
            pop <- 'inuit'
            localGenoSum <- read.table(
                paste('trunc.sim/sim',sim,admix,pop,'genoSum',sep='.')
                ,h=T,as.is=T)
            
            localMerge <- merge(localGenoSum,phe,by=c('FID','IID'))
            tmpLocal <- c(h,var(localGenoSum$genoBetaSum)/var(localMerge$phe),admix,pop,nrow(localMerge),sim)
            names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')
            
            localH <- rbind(localH,tmpLocal)
            
       }
        
        for(pop in c('inuit','danes')){
            for(admix in paste(seq(0.5,0.7,0.1),sep='')){               
                localGenoSum <- read.table(
                    paste('trunc.sim/sim',sim,admix,pop,'genoSum',sep='.')
                    ,h=T,as.is=T)

                localMerge <- merge(localGenoSum,phe,by=c('FID','IID'))
                tmpLocal <- c(h,var(localGenoSum$genoBetaSum)/var(localMerge$phe),admix,pop,nrow(localMerge),sim)
                names(tmpLocal) <- c('simH','localH','admixProp','pop','N','sim')
                
                localH <- rbind(localH,tmpLocal)
                
            }
            
        }
    }
}



t <- data.frame(
    simH=as.numeric(localH[,'simH']),
    localH=as.numeric(localH[,'localH']),
    admixProp=localH[,'admixProp'],
    pop=localH[,'pop'],
    N=as.numeric(localH[,'N']),
    sim=as.numeric(localH[,'sim']),
    stringsAsFactors=F)


Nsum <- sum(t[which(t$pop!='all' & grepl('to',t$admixProp,fixed=T) & t$sim==1 & t$simH==0.4),'N'])

#####
estimates <- NULL
for(h in c(0.4,0.6,0.8)){
    for(sim in 1:100){
        
        if(paste('sim',sim,'h',h,'final.out',sep='.') %in% dir('SNPFAM')){
            snpfam <- read.table(
                paste('SNPFAM/sim',sim,'h',h,'final.out',sep='.'),
                h=F,as.is=T)       
        }else{
            snpfam <- data.frame(V1=NA,V2=NA)
        }

        if(paste('sim',sim,'h',h,'classicPED.out',sep='.') %in% dir('classicPED')){
            ped <- read.table(
                paste('classicPED/sim',sim,'h',h,'classicPED.out',sep='.'),
                h=T,as.is=T)
        }else{
            ped <- data.frame(h2G=NA,Varh2G=NA)
        }

        if(paste('sim',sim,'h',h,'classicPEDwPCA.out',sep='.') %in% dir('classicPEDwPCA')){
            pedwPCA <- read.table(
                paste('classicPEDwPCA/sim',sim,'h',h,'classicPEDwPCA.out',sep='.'),
                h=T,as.is=T)
        }else{
            pedwPCA <- data.frame(h2G=NA,Varh2G=NA)
        }
        
        tmpEst <- c(ped$h2G,sqrt(ped$Varh2G),pedwPCA$h2G,sqrt(pedwPCA$Varh2G),snpfam$V1,snpfam$V2,h,sim)
        names(tmpEst) <- c('ped.h2','ped.SE','pedwPCA.h2','pedwPCA.SE','SNPFAM.h2','SNPFAM.SE','simH','sim')
        estimates <- rbind(estimates,tmpEst)      
    }
}

estDF <- data.frame(
    ped.h2=as.numeric(estimates[,'ped.h2']),
    ped.SE=as.numeric(estimates[,'ped.SE']),
    pedwPCA.h2=as.numeric(estimates[,'pedwPCA.h2']),
    pedwPCA.SE=as.numeric(estimates[,'pedwPCA.SE']),
    SNPFAM.h2=as.numeric(estimates[,'SNPFAM.h2']),
    SNPFAM.SE=as.numeric(estimates[,'SNPFAM.SE']),
    simH=as.numeric(estimates[,'simH']),
    sim=as.numeric(estimates[,'sim']),
    stringsAsFactors=F)




metal <- function(hs,ses){
    w <- 1/(ses^2)
    b <- sum(hs*w)/sum(w)
    se <- sqrt(1/sum(w))
    c(b,se)
}

estCollected <- NULL
for(h in c(0.4,0.6,0.8)){
    #Seperate into hs
    a <- estDF[estDF$simH==h,]
    a <- a[complete.cases(a),]
    
    #Estimates
    pedtmp <- c(metal(a$ped.h2,a$ped.SE),sd(a$ped.h2),'ped',h,NA)
    names(pedtmp) <- c('hLocal','SE','SD','localSub','simH','N')

    pedtmpwPCA <- c(metal(a$pedwPCA.h2,a$pedwPCA.SE),sd(a$pedwPCA.h2),'pedwPCA',h,NA)
    names(pedtmpwPCA) <- c('hLocal','SE','SD','localSub','simH','N')

    snpfamtmp <- c(metal(a$SNPFAM.h2,a$SNPFAM.SE),sd(a$SNPFAM.h2),'snpfam',h,NA)
    names(snpfamtmp) <- c('hLocal','SE','SD','localSub','simH','N')
    
    estCollected <- rbind(estCollected,pedtmp,pedtmpwPCA,snpfamtmp)
 
}


#####
inuitSpecial <- c('0.9.to.0.99','0.99.to.1')
    
locals <- NULL
for(h in c(0.4,0.6,0.8)){
      
    #Seperate into hs
    a <- t[t$simH==h,]
    meanLocalH <- NULL
    #Inuit
    for(i in seq(0.5,0.8,0.1)){
        j <- i+0.1
        sub <- paste(i,'to',j,sep='.')
        pop <- 'inuit'
        N <- a[which(a$admixProp==sub & a$pop==pop),'N']
        cand <- a[which(a$admixProp==sub & a$pop==pop),'localH']

        meanLocalH <- cbind(meanLocalH,cand*(N/Nsum))

        tcand <- t.test(cand)
        
        tmpLocal <- c(tcand$estimate[[1]],
                      (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                      sd(cand),
                      paste('(',i*100,'%-',j*100,'%].',pop,sep=''),
                      h,N[1])
        
        names(tmpLocal) <- c('h2'
                             ,'SE'
                             ,'SD'
                             ,'localSub','simH','N')

        locals <- rbind(locals,tmpLocal)
    }

    #Inuit special
    for(sub in inuitSpecial){
        pop <- 'inuit'
        N <- a[which(a$admixProp==sub & a$pop==pop),'N']
        cand <- a[which(a$admixProp==sub & a$pop==pop),'localH']

        meanLocalH <- cbind(meanLocalH,cand*(N/Nsum))

        tcand <- t.test(cand)
        
        tmpLocal <- c(tcand$estimate[[1]],
                      (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                      sd(cand),
                      ifelse(sub==inuitSpecial[1],
                             '(90%-99%].inuit','(99%-100%].inuit'),
                      h,N[1])
        
        names(tmpLocal) <- c('h2'
                             ,'SE'
                             ,'SD'
                             ,'localSub','simH','N')

        locals <- rbind(locals,tmpLocal)
         }


      #Danes
    for(i in seq(0.5,0.7,0.1)){
        j <- i+0.1
        sub <- paste(i,'to',j,sep='.')
        pop <- 'danes'
        N <- a[which(a$admixProp==sub & a$pop==pop),'N']
        cand <- a[which(a$admixProp==sub & a$pop==pop),'localH']

        meanLocalH <- cbind(meanLocalH,cand*(N/Nsum))
                
        tcand <- t.test(cand)
        
        tmpLocal <- c(tcand$estimate[[1]],
                      (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                      sd(cand),
                      paste('(',i*100,'%-',j*100,'%].',pop,sep=''),h,N[1])
        
        names(tmpLocal) <- c('h2'
                             ,'SE'
                             ,'SD'
                             ,'localSub','simH','N')
       
        locals <- rbind(locals,tmpLocal)
    }

    #Average weigthed local heritabilites
    cand <- apply(meanLocalH,1,sum)
    tcand <- t.test(cand)
        tmpLocal <- c(tcand$estimate[[1]],
                      (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                      sd(cand),
                      'weigthed local avr'
                      ,h,Nsum)
    
    names(tmpLocal) <- c('h2'
                         ,'SE'
                         ,'SD'
                         ,'localSub','simH','N')
    locals <- rbind(locals,tmpLocal)
    
    #All combinded (global h)
    cand <- a[a$pop=='all','localH']
    tcand <- t.test(cand)
    N <- a[which(a$pop=='all'),'N']    
    tmpLocal <- c(tcand$estimate[[1]],
                  (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                  sd(cand),
                  'all'
                  ,h,N[1])
    
    names(tmpLocal) <- c('h2'
                         ,'SE'
                         ,'SD'
                         ,'localSub',
                         'simH','N') 
    
    locals <- rbind(locals,tmpLocal)
    

    ##Thresholds 
    #Special case for super inuit (80%, 90% and 99%, for which there is no danes equilivant. 
    for(sub in c('0.8','0.9','0.99')){
        pop <- 'inuit'       
        cand <- a[which(a$admixProp==sub & a$pop==pop),'localH']
        tcand <- t.test(cand)
        N <- a[which(a$admixProp==sub & a$pop==pop),'N']
        
        tmpLocal <- c(tcand$estimate[[1]],
                      (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                      sd(cand),
                      paste(sub,pop,sep='.'),h,N[1])
        
        names(tmpLocal) <- c('h2'
                             ,'SE'
                             ,'SD'
                             ,'localSub','simH','N')
        
        locals <- rbind(locals,tmpLocal)
        
    }
    
    for(pop in c('inuit','danes')){
        for(sub in paste(seq(0.5,0.7,0.1),sep='')){

            cand <- a[which(a$admixProp==sub & a$pop==pop),'localH']
            tcand <- t.test(cand)
            N <- a[which(a$admixProp==sub & a$pop==pop),'N']
            
            tmpLocal <- c(tcand$estimate[[1]],
                          (tcand$conf.int[[2]]-tcand$estimate[[1]])/1.92,
                          sd(cand),
                          paste(sub,pop,sep='.'),h,N[1])
            
            names(tmpLocal) <- c('h2'
                                 ,'SE'
                                 ,'SD'
                                 ,'localSub','simH','N')
            
            locals <- rbind(locals,tmpLocal)
        }
    }
}


#######
locals <- rbind(locals,estCollected)

localDF <- data.frame(
    hLocal=as.numeric(locals[,'h2']),
    SE=as.numeric(locals[,'SE']),
    SD=as.numeric(locals[,'SD']),
    type=locals[,'localSub'],
    simH=as.numeric(locals[,'simH']),
    N=as.numeric(locals[,'N']),
    stringsAsFactors=T)


localDF$type <- factor(localDF$type,c(
    paste('(',seq(0.5,0.8,0.1)*100,'%-',seq(0.6,0.9,0.1)*100,'%].','inuit',sep=''), #inuit
    '(90%-99%].inuit','(99%-100%].inuit', #inuit specials
    paste('(',seq(0.5,0.7,0.1)*100,'%-',seq(0.6,0.8,0.1)*100,'%].','danes',sep=''), #danes

    paste(seq(0.5,0.9,0.1),'inuit',sep='.'), #threshold inuit
    '0.99.inuit', #threshold inuit specials
    paste(seq(0.5,0.7,0.1),'danes',sep='.'), #threshold danes
    
    'all', 'weigthed local avr', #global and average local sim h
    'ped','pedwPCA','snpfam')) #estimators

localDF$simH <- factor(localDF$simH,c(0.4,0.6,0.8))
localDF$xloc <- (as.integer(localDF$simH)-1)*22+(as.integer(localDF$type)-1)
localDF <- localDF[order(localDF$xloc),]
    
write.table(localDF,'all.different.sim.into.one.table',row.names=F,quote=F)

#####
library(plotrix)
subCol <- c(rep(c('cyan',paste('cyan',1:4,sep=''),'blue3','coral','coral2','coral4'),2),'red','green',1,1,1)

for(h in c(0.4,0.6,0.8)){
    png(paste('comparing.global.with.local.across.100.different.sim.for.h',h,'png',sep='.'),width=800)

    plotDF <- localDF[localDF$simH==h,]
    
    par(mar=c(8,4,4,12)+0.1)
    par(xpd=FALSE)
    
    plotCI(x=plotDF$xloc,y=plotDF$hLocal,
       ui=plotDF$hLocal + 1.92*plotDF$SE,
       li=plotDF$hLocal - 1.92*plotDF$SE,
       pch=c(rep(21,9),rep(22,9),23,23,24,24,25),
       pt.bg=subCol[as.integer(plotDF$type) %% 24],      
       ylab='h2',ylim=c(0,1),main=h,
       xaxt='n',xlab='')


    plotCI(x=plotDF$xloc,y=plotDF$hLocal,
       ui=plotDF$hLocal + 1.92*plotDF$SD,
       li=plotDF$hLocal - 1.92*plotDF$SD,
       pch=c(rep(21,9),rep(22,9),23,23,24,24,25),
       pt.bg=subCol[as.integer(plotDF$type) %% 24],      
       add=T)
    
    axis(side=1,
         at=plotDF$xloc,
         labels=as.character(plotDF$type),
         las=2
         )

    abline(h=h)
    
    par(xpd=TRUE)
    legend(x=max(plotDF$xloc)+1,y=1.075,
           pch=c(rep(21,9),rep(22,9),23,23,24,24,25),
           pt.bg=subCol[as.integer(plotDF$type) %% 24],      
           legend=c(paste(plotDF[1:20,'type'],' (',plotDF[1:20,'N'],')',sep=''),as.character(plotDF[21:23,'type'])),
           bg='white')
    
    dev.off()
}

q()
