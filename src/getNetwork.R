#X: log expression [genes,samples]
# X <- matrix(rnorm(gN*sN,10),nrow=gN,ncol=sN,dimnames=list(paste0("g",1:gN),paste0("s",1:sN)))
getNetwork <- function(X,cor_cutoff=0.7,p_cutoff=0.05,variableN=3000,edge_max=2e6,edge_min=2e3,core=2){
    suppressMessages(require(BiocParallel))
    suppressMessages(require(tidyverse))
    if(nrow(X)>variableN){
        X <- X[order(apply(X,1,sd),decreasing=T)[1:variableN],]
    }
    if(nrow(X)>6e5){## the parallel is much slower than Hmisc
        grpN <- 100
        grp <- c("grp001",rep(paste0("grp",gsub(" ","0",format(1:grpN,width=3))),
                            diff(round(seq(1,nrow(X),length.out=grpN+1)))))
        system.time(res <- list(r=cor(t(X)),
                        P=do.call(rbind,bplapply(split(as.data.frame(X),grp),
                                   geneCorP,X,BPPARAM=MulticoreParam(core,tasks=grpN)))))
    }else{
        system.time(cor_res <- Hmisc.rcorr(as.matrix(t(X))))
    }

    cormat <- cor_res$r
    pmat <- cor_res$P
    ut <- upper.tri(cormat)
    network <- tibble (
        from = rownames(cormat)[row(cormat)[ut]],
        to = rownames(cormat)[col(cormat)[ut]],
        cor  = signif(cormat[ut], 2),
        p = signif(pmat[ut], 2),
        direction = as.integer(sign(cormat[ut]))
    )
    selEdge <- sum(!is.na(network$cor) & abs(network$cor) > cor_cutoff & network$p < p_cutoff)
    #check network size. If it has>5 million rows, use higher cor and lower p to futher reduce size
    if(selEdge>edge_max){
        cor_cutoff <- sort(abs(network$cor),decreasing=T)[edge_max]
        p_cutoff <- sort(network$p)[edge_max]
    }else if(selEdge<edge_min){
        cor_cutoff <- sort(abs(network$cor),decreasing=T)[edge_min]
        p_cutoff <- sort(network$p)[edge_min]
    }
    network <- network %>% mutate_if(is.factor, as.character) %>%
        dplyr::filter(!is.na(cor) & abs(cor) > cor_cutoff & p < p_cutoff)
    return(network)
}

geneCorP <- function(x,y,...){
    x <- as.matrix(x)
    y <- as.matrix(y)
    P <- sapply(split(y,1:nrow(y)),function(yy,X){
        sapply(split(X,1:nrow(X)),function(x1,y1)return(cor.test(x1,y1)$p.value),yy)
    },x)
    dimnames(P) <- list(rownames(x),rownames(y))
    
    return(P)
}
