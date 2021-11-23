
# the input of the expression counts
# 
# input parameters
# X: the count matrix [genes,samples]
# effeL: the effective_length matrix
# batchX: The table of covariates [samples,covariates]
# method: limma, combat_seq(only support one factor variate)
# prior: the small value added to estimated TPM before log2 transform
# minCount, minSample: genes not satisify this will not be adjusted
# 
# output
# the log2 scale of batch corrected TPM
covariateRM <- function(X,effeL,batchX=NULL,method='limma',
                        prior=0.25,minCount=1,minSample=1){
    ## check the input data
    if(sum(!rownames(X)%in%rownames(effeL))>0 || sum(!colnames(X)%in%colnames(effeL))>0)
        stop("sample names and/or gene names does not match between the count table and effective_length table!")
    effeL <- as.matrix(effeL[rownames(X),colnames(X)])
    
    if(!is.null(method) && !is.null(batchX)){
        if(sum(colnames(X)!=rownames(batchX))>0)
            stop("Sample names does not match between the count table and meta information")
        batchX <- batchX[colnames(X),,drop=F]
    }
    # filter low express genes
    X <- as.matrix(X)
    ix <- apply(X,1,function(x)return(sum(x>=minCount)))>=minSample
    sizeF<-NULL
    ## batch removal counts/cpm
    if(is.null(method) || is.null(batchX)){
        ad_X <- X
    }else if(method=="limma"){
        ad_X <- limmaRM(X[ix,],batchX,prior)
        ad_X <- rbind(ad_X,X[!ix,])
        sizeF <- getSizeF(ad_X)
    }else if(method=="combat_seq"){
        ad_X <- ComBatRM(X[ix,],batchX[,1])
        ad_X <- rbind(ad_X,X[!ix,])
        sizeF <- getSizeF(ad_X)
    }else{
        stop("unknown batch removal method!")
    }
    logTPM <- estTPM(ad_X,effeL,sizeF,prior)
    message("Finished TPM estimation!")
    return(logTPM)
}

# suggested from https://support.bioconductor.org/p/121523/
limmaRM <- function(X,batchX,prior=0.25){
    if(nrow(batchX)!=ncol(X))
        stop(paste("Sample number is inconsistent between expression (",ncol(X),") and batch (",nrow(batchX),")"))
    Xmax <- 10^ceiling(log(max(X),10))
    factorX <- c()
    design <- data.frame(matrix(1, nrow(batchX), 1))
    for(i in colnames(batchX)){
        if(is.factor(batchX[,i]) || is.character(batchX[,i])){
            batch <- as.factor(batchX[,i])
            batchX <- batchX[,colnames(batchX)!=i,drop=F]
            contrasts(batch) <- contr.sum(levels(batch))
            batch <- model.matrix(~batch)[, -1, drop = FALSE]
            colnames(batch) <- gsub("batch",i,colnames(batch))
            factorX <- cbind(factorX,batch)
        }
    }
    if(!is.null(factorX)) design <- cbind(design,factorX)
    design <- cbind(design,batchX)

    message("Starting limma batch RM")
    M <- edgeR::DGEList(counts=X)
    M <- edgeR::calcNormFactors(M,method="TMM")
    logCPM <- edgeR::cpm(M,log = TRUE,prior.count = prior)
    ## modified from limma::removeBatchEffect
    fit <- limma::lmFit(logCPM,design)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    logCPM <- as.matrix(logCPM) - beta %*% t(design[,-1,drop=F])
    
    #logCPM <- limma::removeBatchEffect(logCPM, batch = batch)
    nlib <- M$samples$lib.size * M$samples$norm.factors
    prior.scaled <- prior * length(nlib) * nlib / sum(nlib)
    estX <- t((2^(t(logCPM)) * (nlib + 2*prior.scaled) / 1e+06) - prior.scaled)
    estX[estX<0] <- 0
    estX[estX>Xmax] <- Xmax
    message("Finished limma batch RM")
    return(estX)
}
ComBatRM <- function(X,batch){
    ad_X <- sva::ComBat_seq(X,batch)
    return(ad_X)
}
getSizeF <- function(X){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=matrix(as.integer(X),nrow=nrow(X),dimnames = dimnames(X)),
                                              colData=data.frame(row.names=colnames(X)),
                                              design=~1)
    return(DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds)))
}
estTPM <- function(X,L,sizeF=NULL,prior=0.25){
    eTPM <- X/L
    eTPM[!is.finite(eTPM)]<-0
    W <- apply(eTPM,2,sum)
    if(!is.null(sizeF)) W <- sizeF*median(W/sizeF)
    eTPM <- sapply(colnames(X),function(a)return(round(eTPM[,a]/W[a]*10^6,2)))
    return(log2(prior+eTPM))
}






