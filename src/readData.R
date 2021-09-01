
# read dnanexus rsem results, TPM, raw counts and effective length

readData <- function(strF){
    D <- read.table(strF,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    D <- D[!grepl("^ERCC",rownames(D)),]
    # remove the TST-prefix and |postfix
    colnames(D) <- sapply(strsplit(sapply(strsplit(colnames(D),"\\|"),
                                          function(x)return(paste(head(x,-1),collapse="|"))),
                                   "_"),function(x)return(paste(x[-1],collapse="_")))
    return(D)
    
}

readQC <- function(strF){
    qc <- read.table(strF,sep="\t",header=T,as.is=T,comment.char="",row.names=1,check.names=F,quote="")
    attr(qc,'oriNames') <- setNames(gsub("^3","x3",gsub("5","x5",gsub("'","p",gsub(" ","_",gsub("\\%","percentage",colnames(qc)))))),
                                    colnames(qc))
    dimnames(qc) <- list(sapply(strsplit(rownames(qc),"_"),function(x)return(gsub(".genome.sorted$","",paste(x[-1],collapse="_")))),
                         attr(qc,'oriNames'))
    qc <- qc[order(rownames(qc)),,drop=F]
    return(qc)
}

matchQCnames <- function(qc,qNames){
    return(attr(qc,'oriNames')[names(attr(qc,'oriNames'))%in%qNames])
}