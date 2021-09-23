
# read dnanexus rsem results, TPM, raw counts and effective length

readData <- function(strF,sName=NULL,gID=NULL){
    D <- read.table(strF,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    D <- D[!grepl("^ERCC",rownames(D)),]
    # remove the TST-prefix and |postfix
    colnames(D) <- sapply(strsplit(sapply(strsplit(colnames(D),"\\|"),
                                          function(x)return(paste(head(x,-1),collapse="|"))),
                                   "_"),function(x)return(paste(x[-1],collapse="_")))
    ## filtering the data
    if(!is.null(sName)){
        if(sum(!sName%in%colnames(D))>0)
            stop(paste0("Samples (",paste(sName[!sName%in%colnames(D)],collapse=","),
                        ") provided in the sample meta information is not listed in ",strF))
        D <- D[,sName,drop=F]
    }
    if(!is.null(gID)){
        gID <- intersect(gID,rownames(D))
        D <- D[gID,]
        message("\t\tFiltering with ",length(gID)," genes")
    }
    return(D)
    
}

readQC <- function(strF,sName=NULL){
    qc <- read.table(strF,sep="\t",header=T,as.is=T,comment.char="",row.names=1,check.names=F,quote="")
    attr(qc,'oriNames') <- setNames(gsub("^3","x3",gsub("5","x5",gsub("'","p",gsub(" ","_",gsub("\\%","percentage",colnames(qc)))))),
                                    colnames(qc))
    dimnames(qc) <- list(sapply(strsplit(rownames(qc),"_"),function(x)return(gsub(".genome.sorted$","",paste(x[-1],collapse="_")))),
                         attr(qc,'oriNames'))
    if(!is.null(sName)){
        if(sum(!sName%in%rownames(qc))>0)
            stop(paste0("Samples (",paste(sName[!sName%in%rownames(qc)],collapse=","),
                        ") provided in the sample meta information is not listed in ",strF))
        qc <- qc[sName,,drop=F]
    }
    #qc <- qc[order(rownames(qc)),,drop=F]
    
    return(qc)
}

matchQCnames <- function(qc,qNames){
    return(attr(qc,'oriNames')[names(attr(qc,'oriNames'))%in%qNames])
}