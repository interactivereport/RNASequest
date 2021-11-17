#source("utility.R",chdir=T)

.libPaths(c(grep("home",.libPaths(),invert=T,value=T),grep("home",.libPaths(),value=T)))
require(data.table)

## EAinit functions -----
source("gtf.gz_to_gene_info.R")
checkInputDir <- function(strInput,strGenome=NULL){
    strInput <- normalizePath(strInput)
    if(!dir.exists(strInput)){
        stop(paste(strInput, "is not a valid path"))
    }
    #for internal data structure
    return(initInternal(strInput,strGenome))
}
initInternal <- function(strInput,strGannotation){
    strSample <- paste0(strInput,"/samplesheet.tsv")
    pInfo <- NULL
    if(file.exists(strSample)){
        pInfo <- list()
        pInfo[["sInfo"]] <- data.frame(fread(strSample),stringsAsFactors=F)
        pInfo$sInfo <- pInfo$sInfo[,apply(pInfo$sInfo,2,function(x)return(sum(!is.na(x))>0))]
        strConfig <- paste0(strInput,"/config.json")
        if(!file.exists(strConfig)) stop("Incomplete DNAnexus download: missing ",strConfig)
        pInfo[["prjID"]] <- pInfo[["prjTitle"]] <- getProjectID(strConfig)
        pInfo[["species"]] <- getSpecies(strConfig)
        pInfo[["gAnnotion"]] <- extractAnnotation(strConfig,strGannotation)
        
        strPrj <- gsub("tsv$","json",strSample)
        if(file.exists(strPrj)){
            prjInfo <- rjson::fromJSON(file=strPrj)$Project
            pInfo[["prjID"]] <- prjInfo$TSTID
            pInfo[["prjTitle"]] <- prjInfo$Study_Title
        }
        
        pInfo[["strCount"]] <- normalizePath(paste0(strInput,"/combine_rsem_outputs/genes.estcount_table.txt"))
        pInfo[["strEffLength"]] <- normalizePath(extractEffectLength(strInput))
        pInfo[["strSeqQC"]] <- normalizePath(paste0(strInput,"/combine_rnaseqc/combined.metrics.tsv"))
    }
    return(pInfo)
}
extractEffectLength<-function(strPath){
    message("Extracting effective length ...")
    strF <- paste0(strPath,"/combine_rsem_outputs/genes.effective_length.txt")
    if(file.exists(strF)) return(strF)
    if(dir.exists(paste0(strPath,"/rsem"))){
        D <- NULL
        selColumn <- c("gene_id","effective_length")
        for(i in list.files(paste0(strPath,"/rsem"),"genes.results.gz",full.names=T)){
            X <- fread(i,sep="\t",header=T)[,..selColumn]#,as.is=T,row.names=1
            colnames(X) <- c(selColumn[1],gsub(".genes.results.gz",paste0("|",selColumn[2]),basename(i)))
            if(is.null(D)) D <- X
            else D <- merge(D,X,by=selColumn[1],all=T)
        }
        write.table(D,file=strF,sep="\t",row.names=F,quote=F)
    }else{
        stop("Cannot create Effective length file, check DNAnexus 'rsem' folder")
    }
    return(strF)
}
extractAnnotation <- function(strF,strPath=NULL){
    message("Create gene annotation ...")
    if(grepl("json$",strF)){
        gConfig <- rjson::fromJSON(file=strF)
        gtfPath <- list.files(paste0(strPath,"/rnaseq/",
                                     gConfig$global_params$reference$species,
                                     "/",gConfig$global_params$reference$version),
                              "gtf.gz$",full.names=T)[1]
        strF <- gsub("gtf.gz$","gene_info.csv",gtfPath)
        if(!file.exists(strF)){
            gtf.gz_to_gene_info(gtfPath)
        }
    }
    #columns in strF: unique ID, gene_name, gene_type,....
    gInfo <- read.csv(strF,as.is=T,check.names=F)
    gInfoName <- colnames(gInfo)
    gInfo <- cbind(0:(nrow(gInfo)-1),gInfo)
    dimnames(gInfo) <- list(gInfo[,2],
                            c("id","UniqueID","Gene.Name","Biotype",tail(gInfoName,-3)))
    ## if the lookup table is available
    strLookup <- paste0(strPath,"/rnaseq/lookup/",
                        gConfig$global_params$reference$version,
                        ".csv")
    if(file.exists(strLookup)){
        gLookup <- read.csv(strLookup,row.names=1,as.is=T,check.names = F)
        gInfo <- merge(gInfo,gLookup,by="row.names",all=T,sort=F)
        rownames(gInfo) <- gInfo[,1]
        gInfo <- gInfo[,-1]
    }
    return(gInfo)
}
getSpecies <- function(strF){
    gConfig <- rjson::fromJSON(file=strF)
    return(gConfig$global_params$reference$species)
}
getProjectID <- function(strF){
    gConfig <- rjson::fromJSON(file=strF)
    pID <- gConfig$global_params$internal_project_id
    if(!is.null(pID)) return(pID)
    for(one in gConfig$execution$manifest){
        pID <- c(pID,one$project)
    }
    return(paste(unique(pID),collapse="_"))
}

createEmptyComp <- function(strComp){
    message("Create empty comparison template ...")
    comTitle <- c("CompareName",
                  "Subsetting_group",
                  "Model",
                  "Covariate_levels",
                  "Group_name",
                  "Group_test",
                  "Group_ctrl",
                  "Analysis_method",
                  "Shrink_logFC",
                  "LFC_cutoff")
    cat(paste(comTitle,collapse=","),"\n",sep="",file=strComp)
}
saveGeneAnnotation <- function(gAnno,strF){
    if(!is.null(gAnno)){
        write.csv(gAnno,file=strF)
    }
}
cleanTST <- function(strSrc,strDest=NULL,sep="\t"){
    if(is.null(strSrc) || !file.exists(strSrc)) return(NULL)
    D <- read.table(strSrc,sep=sep,header=T,as.is=T,check.names=F,quote="")
    if(sum(grepl("^TST",colnames(D)))>0){
        colnames(D) <- sapply(strsplit(sapply(strsplit(colnames(D),"\\|"),head,1),
                                       "_"),
                              function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
    }else if(sum(grepl("^TST",D[,1]))>0){
        D[,1] <- sapply(strsplit(gsub(".genome.sorted","",D[,1]),
                                 "_"),
                        function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
    }
    if(!is.null(strDest)) write.table(D,strDest,row.names=F,sep=sep)
    else return(D)
}
appendMeta <- function(pInfo,sample_name,selQC){
    message("Appending sequencing QC into sample meta file ...")
    if(!is.null(pInfo$sInfo) && !is.null(pInfo$strSeqQC)){
        qc <- cleanTST(pInfo$strSeqQC)
        rownames(qc) <- qc[,1]
        qc <- qc[,selQC,drop=F]
        colnames(qc) <- gsub("^3","x3",gsub("5","x5",gsub("'","p",gsub(" ","_",gsub("\\%","percentage",colnames(qc))))))
        rownames(pInfo$sInfo) <- pInfo$sInfo[,sample_name]
        if(sum(!rownames(pInfo$sInfo)%in%rownames(qc))>0)
            stop(paste0("Sample Names (",
                       paste(rownames(pInfo$sInfo)[!rownames(pInfo$sInfo)%in%rownames(qc)],collapse=", "),
                       ") specified in sample sheet cannot be found in seqQC file ",
                       pInfo$strSeqQC))
        meta <- merge(pInfo$sInfo,qc,by="row.names",sort=F)
        rownames(meta) <- meta[,1]
        meta <- meta[,-1]
        pInfo$sInfo <- meta
    }
    return(pInfo)

}
getCovariates <- function(pInfo,notCovariates=NULL){
    if(!is.null(pInfo$sInfo)){
        covariates <- NULL
        for(i in setdiff(colnames(pInfo$sInfo),notCovariates)){
            if(length(unique(pInfo$sInfo[,i]))>1) covariates <- c(covariates,i)
        }
        pInfo$covariates <- covariates
    }
    return(pInfo)
    
}
createInit <- function(strInput,configTmp,pInfo){
    strInput <- normalizePath(strInput)
    message("Creating project folder")
    ix <- 0
    while(dir.exists(strOut<-paste0(strInput,"/EA",gsub("\\-","",Sys.Date()),"_",ix))){
        ix <- ix+1
    }
    system(paste("mkdir -p",paste0(strOut,"/data")))
    
    strCount <- paste0(strOut,"/data/count.tsv")
    strEffLength <- paste0(strOut,"/data/effLength.tsv")
    strSeqQC <- paste0(strOut,"/data/seqQC.tsv")
    
    strMeta <- paste0(strOut,"/data/sampleMeta.csv")
    strMetaFactor <- paste0(strOut,"/data/sampleMetaFactor.yml")
    strGinfo <- paste0(strOut,"/data/geneAnnotation.csv")
    saveGeneAnnotation(pInfo$gAnnotion,strGinfo)
    
    strComp <- paste0(strOut,"/data/compareInfo.csv")
    createEmptyComp(strComp)
    if(is.null(pInfo)){
        configTmp <- gsub("initPrjName"," #required",configTmp)
        configTmp <- gsub("initPrjTitle"," #optional",configTmp)
        configTmp <- gsub("prj_counts"," #required",configTmp)
        configTmp <- gsub("prj_effLength"," #if provided prj_TPM is ignored",configTmp)
        configTmp <- gsub("prj_seqQC"," #optional",configTmp)
        configTmp <- gsub("prj_TPM","  #prj_effLength or prj_TPM is required",configTmp)
        configTmp <- gsub("initPrjMeta"," #required",configTmp)
        configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
        configTmp <- gsub("initSpecies"," #required",configTmp)
        configTmp <- gsub("initGeneAnnotation"," #optional",configTmp)
        configTmp <- gsub("initOutput",strOut,configTmp)
        configTmp <- gsub("initCovariates","",configTmp)
        configTmp <- gsub("initPrjComp",strComp,configTmp)
    }else{
        cleanTST(pInfo[["strCount"]],strDest=strCount)
        cleanTST(pInfo[["strEffLength"]],strDest=strEffLength)
        cleanTST(pInfo[["strSeqQC"]],strDest=strSeqQC)
        write.csv(pInfo$sInfo,file=strMeta,row.names=F)
        
        configTmp <- gsub("initPrjName",pInfo[["prjID"]],configTmp)
        configTmp <- gsub("initPrjTitle",pInfo[["prjTitle"]],configTmp)
        configTmp <- gsub("initCounts",strCount,configTmp)
        configTmp <- gsub("initEffLength",strEffLength,configTmp)
        configTmp <- gsub("initSeqQC",strSeqQC,configTmp)
        configTmp <- gsub("initTPM","",configTmp)
        
        configTmp <- gsub("initPrjMeta",strMeta,configTmp)
        configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
        configTmp <- gsub("initSpecies",pInfo$species,configTmp)
        configTmp <- gsub("initGeneAnnotation",strGinfo,configTmp)
        configTmp <- gsub("initOutput",strOut,configTmp)
        configTmp <- gsub("initCovariates",paste0("[",paste(pInfo$covariates,collapse=","),"]"),configTmp)
        configTmp <- gsub("initPrjComp",strComp,configTmp)
        
        configTmp <- gsub("^shinyOne_Title:",paste("shinyOne_Title:",pInfo$prjID),configTmp)
        configTmp <- gsub("^shinyOne_Description:",paste0("shinyOne_Description:",pInfo$prjTitle),configTmp)
        
    }
    cat(paste(configTmp,collapse="\n"),"\n",sep="",file=paste0(strOut,"/config.yml"))
    return(list(strOut=strOut,strComp=strComp))
}
finishInit <- function(strMsg){
    message("==========================================")
    message("ExpressionAnalysis project folder is created at ",strMsg$strOut)
    message("-----> 'EAqc' can be used to identify the covariates to be adjusted as:")
    message("\t\tEAqc ",strMsg$strOut,"/config.yml\n\n")
    message("----->'EArun' can be used to obtain the QuickOmics objects after comparison definition file is updated:")
    message("\t\t\t",strMsg$strComp)
    message("\t\tEArun ",strMsg$strOut,"/config.yml\n\n")
    message("-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.\n")
    
    message("Powered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")
}

## Reading/Checking ALL EA input data -----
checkConfig <- function(config){
    if(is.null(config$sample_name)) stop("sample_name is required in the config. Default is Sample_Name")
    if(config$prj_name!="initPrjName"){
        if(is.null(config$prj_counts) || !file.exists(config$prj_counts))
            stop(paste("The count file is required, ",config$prj_counts,", does NOT exist!"))
        if((is.null(config$prj_effLength) || !file.exists(config$prj_effLength)) && 
           (is.null(config$prj_TPM) || !file.exists(config$prj_TPM)))
            stop(paste0("The effective length (",config$prj_effLength,") or TPM (",config$prj_TPM,") is required!"))
        if(!is.null(config$prj_seqQC) && !file.exists(config$prj_seqQC))
            stop(paste("The prj_seqQC (",config$prj_seqQC,") does NOT exist!"))
        
        if(!file.exists(config$sample_meta)) stop(paste("The sample meta file, ",config$sample_meta,", does NOT exist!"))
        if(is.null(config$species)) stop("species is required")
        if(!is.null(config$gene_annotation) && !file.exists(config$gene_annotation))
            stop(paste("The gene annotation file, ",config$gene_annotation,", does NOT exist!"))
        if(!file.exists(config$comparison_file))
            stop(paste("The comparison definition file, ",config$comparison_file,", does NOT exist!"))
    }
}
getEAData <- function(config){
    D <- getMeta(config)
    D <- c(D,getCounts(config,rownames(D$meta)))
    D <- c(D,getEffLength(config,colnames(D$counts),rownames(D$counts)))
    D <- c(D,getTPM(config,colnames(D$counts),rownames(D$counts)))
    D <- c(D,getSeqQC(config,colnames(D$counts)))
    if(is.null(D$logTPM)){
        D$logTPM <- covariateRM(D$counts,D$effLength,prior=config$count_prior)
    }
    if(is.null(D$gInfo)){
        D$gInfo <- data.frame(row.names = rownames(D),
                              id=0:(nrow(D)-1),
                              UniqueID=rownames(D),
                              Gene.Name=rownames(D))
    }
    return(D)
}
getMeta <- function(config){
    message("reading sample meta")
    meta <- read.csv(config$sample_meta,check.names=F,as.is=T)
    checkMeta(meta,config)
    rownames(meta) <- meta[,config$sample_name]
    return(list(meta=meta))
}
checkMeta <- function(meta,config){
    message("checking against config file")
    if(!config$sample_name%in%colnames(meta))
        stop(paste0("sample_name (",config$sample_name,
                    ") is NOT a column in sample meta file (",
                    config$sample_meta,")!"))
    if(sum(duplicated(meta[,config$sample_name])))
        stop(paste0("sample_name column (",config$sample_name,") contains duplicates in sample meta file."))
    if(!is.null(config$sample_alias)){
        if(config$sample_alias%in%colnames(meta))
            stop(paste0("sample_alias (",config$sample_alias,") is NOT a column in the sample meta file"))
        if(sum(duplicated(meta[,config$sample_alias]))>0)
            stop(paste0("sample_alias column (",config$sample_alias,") contains duplicates in the sample meta file"))
    }
    if(!is.null(config$split_meta) && sum(config$split_meta==colnames(meta))!=1)
        stop(paste0("split_meta (",config$split_meta,") in config is NOT defined in the sample meta file."))
    
    if(!is.null(config$covariates_check) && sum(!config$covariates_check%in%colnames(meta))>0)
        stop(paste0("covariates_check variables defined in config (",
                    paste(config$covariates_check[!config$covariates_check%in%colnames(meta)],collapse=", "),
                    "), is NOT defined in the sample meta file"))
    
    if(!is.null(config$covariates_adjust)){
        if( sum(!config$covariates_adjust%in%colnames(meta))>0)
            stop(paste0("covariates_adjust variables defined in config (",
                        config$covariates_adjust[!config$covariates_adjust%in%colnames(meta)],
                        "), is NOT defined in the sample meta file"))
        for(one in config$covariates_adjust){
            if(length(unique(meta[,one]))<2)
                stop(paste0("covariates_adjust (",one,") only contains one unique value in the sample meta and cannot be used for adjusting"))
        }
    }
}
getCounts <- function(config,sID){
    message("reading sample counts")
    D <- read.table(config$prj_counts,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    ix <- apply(as.matrix(D),1,function(x)return(sum(x>=config$min_count)))>=config$min_sample
    message("\tFiltering genes (",sum(ix),") with minimal counts ",
            config$min_count," in at least ",config$min_sample," samples")
    D <- D[ix,]
    gInfo <- NULL
    if(!is.null(config$gene_annotation)){
        gInfo <-  read.csv(config$gene_annotation,row.names=1,as.is=T)
        gID <- intersect(rownames(D),rownames(gInfo))
        if(length(gID)<nrow(D)){
            D <- D[gID,]
            gInfo <- gInfo[gID,]
            message("\tFiltering genes (",length(gID),") by gene annotation")
        }
    }
    D <- checkSampleName(D,sID)
    return(list(counts=D,gInfo=gInfo))
}
checkSampleName <- function(D,sID){
    if(sum(!sID%in%colnames(D))>0)
        stop(paste0("samples (",paste(sID[!sID%in%colnames(D)],collapse=","),
                    ") defined in sample meta table are NOT available in count matrix"))
    D <- D[,sID,drop=F]
    return(D)
}
checkGeneName <- function(D,gID){
    if(sum(!gID%in%rownames(D))>0)
        stop(paste0("genes defined in count table (",length(gID),") are NOT available in effective length or TPM table (",nrow(D),")"))
    D <- D[gID,,drop=F]
    return(D)
}
getEffLength <- function(config,sID,gID){
    D <- NULL
    if(!is.null(config$prj_effLength) && file.exists(config$prj_effLength)){
        message("reading effective length")
        D <- read.table(config$prj_effLength,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
        D <- checkSampleName(D,sID)
        D <- checkGeneName(D,gID)
    }
    return(list(effLength=D))
}
getTPM <- function(config,sID,gID){
    D <- NULL
    if(!is.null(config$prj_TPM) && file.exists(config$prj_TPM)){
        message("reading sample TPM")
        D <- read.table(config$prj_TPM,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
        D <- checkSampleName(D,sID)
        D <- checkGeneName(D,gID)
        D <- log2(config$count_prior+D)
    }
    return(list(logTPM=D))
}
getSeqQC <- function(config,sID){
    if(is.null(config$prj_seqQC)) return(NULL)
    message("reading sequence QC")
    D <- read.table(config$prj_seqQC,sep="\t",header=T,as.is=T,check.names=F,row.names=1)
    if(sum(!sID%in%rownames(D))>0)
        stop(paste0("samples (",paste(sID[!sID%in%colnames(D)],collapse=","),
                    ") defined in sample meta table are NOT available in sequence QC table"))
    D <- D[sID,,drop=F]
    return(list(seqQC=D))
}
useAlias <- function(config,D){
    if(!is.null(config$sample_alias)){
        message("Applying alias")
        colnames(D$counts) <- rownames(D$meta) <- D$meta[,config$sample_alias]
        if(!is.null(D$effLength)) colnames(D$effLength) <- rownames(D$meta)
        if(!is.null(D$TPM)) colnames(D$TPM) <- rownames(D$meta)
        if(!is.null(D$seqQC)) colnames(D$seqQC) <- rownames(D$meta)
    }
    return(D)
}

## batch removal based on counts and effective length --------
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
        ad_X <- covariateRM_limmaRM(X[ix,],batchX,prior)
        ad_X <- rbind(ad_X,X[!ix,])
        sizeF <- covariateRM_getSizeF(ad_X)
    }else if(method=="combat_seq"){
        ad_X <- covariateRM_ComBatRM(X[ix,],batchX[,1])
        ad_X <- rbind(ad_X,X[!ix,])
        sizeF <- covariateRM_getSizeF(ad_X)
    }else{
        stop("unknown batch removal method!")
    }
    logTPM <- covariateRM_estTPM(ad_X,effeL,sizeF,prior)
    message("Finished TPM estimation!")
    return(logTPM)
}

# suggested from https://support.bioconductor.org/p/121523/
covariateRM_limmaRM <- function(X,batchX,prior=0.25){
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
covariateRM_ComBatRM <- function(X,batch){
    ad_X <- sva::ComBat_seq(X,batch)
    return(ad_X)
}
covariateRM_getSizeF <- function(X){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=matrix(as.integer(X),nrow=nrow(X),dimnames = dimnames(X)),
                                          colData=data.frame(row.names=colnames(X)),
                                          design=~1)
    return(DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds)))
}
covariateRM_estTPM <- function(X,L,sizeF=NULL,prior=0.25){
    eTPM <- X/L
    eTPM[!is.finite(eTPM)]<-0
    W <- apply(eTPM,2,sum)
    if(!is.null(sizeF)) W <- sizeF*median(W/sizeF)
    eTPM <- sapply(colnames(X),function(a)return(round(eTPM[,a]/W[a]*10^6,2)))
    return(log2(prior+eTPM))
}

## EAqc functions ------
require(ggplot2)
require(reshape2)
plotAlignQC <- function(estT,strPDF,estC=NULL,qc=NULL,prioQC=NULL,topN=c(1,10,30),gInfo=NULL,replot=F){#,50,100
    if(file.exists(strPDF) && !replot) return()
    message("plotting sequencing QC @",strPDF)
    if(!is.null(gInfo) & sum(rownames(gInfo)!=gInfo$Gene.Name)>0){
        rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
        if(!is.null(estC)) rownames(estC) <- paste(rownames(estC),gInfo[rownames(estC),"Gene.Name"],sep="|")
    }
    
    pdfW <- max(ncol(estT)/10+2,6)
    pdf(strPDF,width=pdfW,height=6)
    ## top genes ratio----
    topN <- setNames(topN,paste0("Top",topN))
    D <- t(apply(estT,2,function(x){
        x <- sort(x,decreasing=T)
        return(sapply(topN,function(i)return(sum(x[1:i])/sum(x)*100)))
    }))
    D <- cbind(sID=factor(rownames(D),levels=rownames(D)),data.frame(D))
    for(i in colnames(D)){
        if(i=="sID") next
        print(ggplot(D,aes_string(x="sID",y=i))+
                  geom_bar(stat="identity")+
                  xlab("")+ylab("% of Total TPM")+
                  ggtitle(paste(i,"genes"))+
                  theme_minimal()+
                  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                        legend.position = "none"))
    }
    ## union top genes across samples ------
    topUnion <- 30
    p <- plotTopGeneRatio(estT,30)
    write.csv(p$data,file=gsub("pdf","unionTop.TPM.csv",strPDF),row.names=F)
    if(pdfW>8) print(p+theme(aspect.ratio=0.75))
    else print(p)
    # same genes for counts
    if(!is.null(estC)){
        pC <- plotTopGeneRatio(estC,30,selG=levels(p$data[,2]))+
            xlab("% of total counts")+ggtitle(p$labels$title)
        if(pdfW>8) print(pC+theme(aspect.ratio=0.75))
        else print(pC)
    }
    ## sequence qc ----
    if(!is.null(qc)){
        ## intergenic, intronic and exonic -----
        selN <- c("Exonic Rate","Intronic Rate","Intergenic Rate")
        if(sum(selN%in%colnames(qc))==length(selN)){
            D = melt(as.matrix(qc[,colnames(qc)%in%selN]))
            D$Var1 <- factor(D$Var1,levels=rownames(qc))
            print(ggplot(D,aes(x=Var1,y=value,fill=Var2))+
                      geom_bar(position="stack",stat="identity")+
                      ylab("Fraction of Reads")+xlab("")+
                      ggtitle("Mapped reads allocation")+
                      ylim(0,1)+theme_minimal()+
                      scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
                      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                            legend.position = "top")+
                      guides(fill=guide_legend(title="")))
            qc <- qc[,!colnames(qc)%in%selN,drop=F]
        }
        ## all rest qc -----
        qc <- cbind(sID=factor(rownames(qc),levels=rownames(qc)),qc)
        selQC <- colnames(qc)%in%prioQC
        colnames(qc) <- make.names(colnames(qc))
        for(i in c(colnames(qc)[selQC],colnames(qc)[!selQC])){
            if(!is.numeric(qc[1,i])) next
            print(ggplot(qc,aes_string(x="sID",y=i))+
                      geom_bar(stat="identity")+
                      xlab("")+ylab("")+
                      ggtitle(i)+
                      theme_minimal()+
                      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                            legend.position = "none"))
        }
    }
    ## -----
    a <- dev.off()
}
plotTopGeneRatio <- function(X,topN,maxN=90,selG=NULL){
    selLevel <- selG
    if(is.null(selG)){
        selG <- unique(as.vector(apply(X,2,function(x)return(names(sort(x,decreasing=T)[1:topN])))))
        while(length(selG)>maxN){
            topN <- topN-1
            selG <- unique(as.vector(apply(X,2,function(x)return(names(sort(x,decreasing=T)[1:topN])))))
        }
    }
    Xsum <- apply(X,2,sum)
    D <- apply(X[selG,],1,function(x)return(100*x/Xsum))
    D <- reshape2::melt(D[,order(apply(D,2,median))])
    if(!is.null(selLevel)) D$Var2 <- factor(D$Var2,selLevel)
    p <- ggplot(D,aes(x=value,y=Var2))+
        geom_point(color="grey50",alpha=0.4,size=1)+
        geom_boxplot(color="#ff7f00",outlier.shape = NA,alpha=0)+
        xlab("% of total TPM")+ylab("")+
        ggtitle(paste("Union of top",topN,"expressed genes"))+
        theme_minimal()+
        theme(axis.text.y=element_text(size=12-(length(selG)-10)/10))
    return(p)
}

plotGeneLength <- function(config,estC,effL=NULL,logTPM=NULL,gInfo=NULL){
    if(!config$geneLength_QC) return()
    if(is.null(effL) && (is.null(gInfo) || !"Length"%in%colnames(gInfo))) return()
    pdf(paste0(config$output,"/geneLengthQC.pdf"))
    par(mar=c(3,3,2,0.5)+0.1,mgp=c(1.1,0.2,0),tcl=-0.1)
    if(!is.null(effL)){
        plotGeneLengthOne(log2(1+cbind(Length=apply(effL[rownames(estC),],1,median),estC)),
                          list(mean=mean,median=median,SD=sd),
                          xlab="log2(median gene effective length +1)",
                          main="est counts")
        if(!is.null(logTPM))
            plotGeneLengthOne(cbind(Length=log2(1+apply(effL[rownames(logTPM),],1,median)),logTPM),
                              list(mean=mean,median=median,SD=sd),
                              xlab="log2(median gene effective length +1)",
                              main="TPM")
    }
    if(!is.null(gInfo) && "Length"%in%colnames(gInfo)){
        plotGeneLengthOne(log2(1+cbind(Length=gInfo[,"Length"],estC)),
                          list(mean=mean,median=median,SD=sd),
                          xlab="log2(GTF gene length +1)",
                          main="est counts")
        if(!is.null(logTPM))
            plotGeneLengthOne(cbind(Length=log2(1+gInfo[,"Length"]),logTPM),
                              list(mean=mean,median=median,SD=sd),
                              xlab="log2(GTF gene length +1)",
                              main="TPM")
    }
    a <- dev.off()

}
plotGeneLengthOne <- function(X,funs,...){
    imageCOL <- c("#FFFFFFFF","#3300FF","#2D1CFF","#2839FF","#2255FF","#1C71FF","#178EFF","#11AAFF",
                  "#0BC6FF","#06E3FF","#00FFFF","#00FFFF","#17FFE3","#2DFFC6","#44FFAA","#5BFF8E",
                  "#71FF71","#88FF55","#9FFF39","#B5FF1C","#CCFF00",
                  "#CCFF00","#D2F400","#D7E800","#DDDD00","#E3D200","#E8C600","#EEBB00","#F4B000",
                  "#F9A400","#FF9900","#FF9900","#F68800","#EC7700","#E36600","#D95500","#D04400",
                  "#C63300","#BD2200","#B31100","#AA0000")
    for(i in names(funs)){
        y <- apply(X[,-1],1,funs[[i]])
        x <- X[,"Length"]
        #x <- X[y>0,"Length"]
        #y <- y[y>0]
        plot(c(),c(),xlim=range(x),ylim=range(y),
             ylab=paste(i,"of log2(Expression+1)"),...)
        index <- rep(T,length(x))
        tryM <- try(f1 <- MASS::kde2d(x,y,n=c(90,60)),silent=T)
        if(!is.null(names(tryM))){
            image(f1,col=imageCOL,add=T)
            imageZero <- diff(range(f1$z))/(length(imageCOL)-1) # ~2% dots
            index <- apply(cbind(x,y),1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
        }
        points(x[index],y[index],pch=20,col=imageCOL[2],cex=1)
    }
}





## general functions used by all ------
saveSeesionInfo <- function(strF,strSRC){
    sink(strF)
    cat("EA version: git-",
        system(paste('git --git-dir',normalizePath(paste0(strSRC,"../.git")),'rev-parse HEAD'),
               intern=T),sep="","\n\n")
    print(sessionInfo())
    sink()
    
}



