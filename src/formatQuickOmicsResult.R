
formatQuickOmicsResult <- function(DEGs,logTPM,grp,gInfo){
    Dw <- data.frame(gInfo[rownames(logTPM),c("UniqueID","Gene.Name",'id')],
                    Intensity=apply(logTPM,1,mean))
    for(i in unique(grp)){
        if(sum(grp==i)<1){
            next
            message("===== warning: no sample for ",i)
        }else if(sum(grp==i)<1){
            tmp <- cbind(logTPM[,grp==i],rep(0,nrow(logTPM)))
            colnames(tmp) <- paste(i,c("Mean","sd"),sep="_")
            Dw <- cbind(Dw,tmp)
            message("===== warning: one sample for ",i)
        }else{
            Dw <- cbind(Dw,t(apply(logTPM[,grp==i,drop=F],1,
                                   function(x)return(setNames(c(mean(x),sd(x)),paste(i,c("Mean","sd"),sep="_"))))) )
        }
    }
    for(i in names(DEGs)){
        Dw <- merge(Dw,DEGs[[i]]$DEG,by="row.names",all=T,sort=F)
        rownames(Dw) <- Dw[,1]
        Dw <- Dw[,-1]
    }
    for(i in grep("logFC$",colnames(Dw))){
        Dw[is.na(Dw[,i]),i] <- 0
    }
    for(i in grep("P.value$",colnames(Dw))){
        Dw[is.na(Dw[,i]),i] <- 1
    }

    Dl <- NULL
    for(i in names(DEGs)){
        res <- DEGs[[i]]$DEG
        colnames(res) <- c("logFC","P.Value","Adj.P.Value")
        res <- cbind(UniqueID=rownames(res),
                     test=factor(i,levels=names(DEGs)),
                     res)
        if(is.null(Dl)){
            Dl <- res
        }else{
            Dl <- rbind(Dl,res)
        }
    }
    return(list(Dw=Dw,Dl=Dl))
}


formatQuickOmicsMeta <- function(meta,comNames){
    MetaData <- list(sampleid=rownames(meta),
                     group=meta$group,
                     Order=unique(meta$group),
                     ComparePairs=comNames)
    MetaData <- cbind(as.data.frame(lapply(MetaData,'length<-',max(sapply(MetaData,length))),stringsAsFactors=F),
                      meta[,-1,drop=F])
    suppressWarnings(MetaData[is.na(MetaData)] <- "")
    return(MetaData)
}