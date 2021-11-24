#source("utility.R",chdir=T)
.libPaths(c(grep("home",.libPaths(),invert=T,value=T),grep("home",.libPaths(),value=T)))

## save session ------
saveSessionInfo <- function(strF,strSRC){
    message("\nPowered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")
    sink(strF)
    cat("EA version: git-",
        system(paste('git --git-dir',normalizePath(paste0(strSRC,"../.git")),'rev-parse HEAD'),
               intern=T),sep="","\n\n")
    print(sessionInfo())
    sink()
    
}
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
    message("saving initialization ...")
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
        configTmp <- gsub("^shinyOne_Description:",paste("shinyOne_Description:",pInfo$prjTitle),configTmp)
        fastrUN <- "FASTR_USER_FULLNAME: "
        configTmp <- gsub("^shinyOne_Data_Generated_By:",
                          paste("shinyOne_Data_Generated_By:",
                                gsub(fastrUN,"",
                                     system(paste0("grep ",fastrUN," ",strInput,"/analysis-*"),
                                            intern=T))),
                          configTmp)
        
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
}

## Reading/Checking ALL EA input data -----
require(data.table)
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
getEAData <- function(config,withCom=F){
    D <- getMeta(config)
    if(withCom)
        D$comp_info <- checkComparisonInfo(read_file(config$comparison_file,T),
                                           D$meta,config$comparison_file)
    D <- c(D,getCounts(config,rownames(D$meta)))
    D <- c(D,getEffLength(config,colnames(D$counts),rownames(D$counts)))
    D <- c(D,getTPM(config,colnames(D$counts),rownames(D$counts)))
    D <- c(D,getSeqQC(config,colnames(D$counts)))
    if(is.null(D$logTPM)){
        D$logTPM <- as.data.frame(covariateRM(D$counts,D$effLength,prior=config$count_prior))
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
    meta <- metaFactor(meta,config$sample_factor)
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
    
    if(!is.null(config$sample_group) && length(config$sample_group)>0){
        if(sum(!config$sample_group%in%colnames(meta))>0){
            stop(paste0("sample_group (",
                        paste(config$sample_group[!config$sample_group%in%colnames(meta)],collapse=","),
                        ") defined in config is NOT included in the sample meta file"))
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
    #print(head(design))
    message("Starting limma batch RM")
    M <- edgeR::DGEList(counts=X)
    M <- edgeR::calcNormFactors(M,method="TMM")
    logCPM <- edgeR::cpm(M,log = TRUE,prior.count = prior)
    ## modified from limma::removeBatchEffect
    fit <- limma::lmFit(logCPM,design,method="robust")
    beta <- fit$coefficients[, -1, drop = FALSE]
    #print(apply(beta,2,summary))
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

## meta factor ----
metaFactor <- function(meta,strMetaFactor,addFactor=NULL){
    if(is.null(strMetaFactor))return(meta)
    ## retrieve the meta factor or initialize one ----
    meta <- metaFactor_checkNAempty(meta)
    if(file.exists(strMetaFactor)){
        metaFactor <- yaml::read_yaml(strMetaFactor)
    }else{
        metaFactor <- metaFactor_addFactor(meta,strMetaFactor)
    }
    ## if any annotation removed from meta but existing in the meta factor ----
    if(sum(!names(metaFactor)%in%colnames(meta))>0){
        warning("The following removed from sampleMeta table, updating metaFactor file:")
        message(paste(names(metaFactor)[!names(metaFactor)%in%colnames(meta)],collapse=","))
        metaFactor <- metaFactor[names(metaFactor)%in%colnames(meta)]
        metaFactor_saveYaml(metaFactor,strMetaFactor)
    }
    ## add additional (new) meta annotation into the meta factor -------
    # since there are numerical annotation, the following mostly always will be executed
    if(sum(!colnames(meta)%in%names(metaFactor))>0)
        metaFactor <- metaFactor_addFactor(meta,strMetaFactor,metaFactor)
    
    ## apply the meta factor to the meta table -----
    for(i in names(metaFactor)){
        ## add additional (new) entry into the meta factor
        oneMetaUnique <- unique(meta[,i])
        if(sum(!oneMetaUnique%in%metaFactor[[i]])>0){
            warning(paste0(paste(oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]],collapse=","),
                           " from ",i," are not defined in meta factor file"))
            message("Please update metaFactor file to avoid this warning message")
            metaFactor[[i]] <- c(metaFactor[[i]],oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]])
        }
        if(sum(!metaFactor[[i]]%in%oneMetaUnique)>0){
            warning(paste0(paste(metaFactor[[i]][!metaFactor[[i]]%in%oneMetaUnique],collapse=","),
                           " from ",i," are missing from sample file"))
            message("Please update metaFactor file to avoid this warning message")
            metaFactor[[i]] <- metaFactor[[i]][metaFactor[[i]]%in%oneMetaUnique]
        }
        if(sum(!oneMetaUnique%in%metaFactor[[i]])>0){
            warning(paste0(paste(oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]],collapse=","),
                           " from ",i," are not defined in meta factor file"))
            message("Please update metaFactor file to avoid this warning message")
            metaFactor[[i]] <- c(metaFactor[[i]],oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]])
        }
        meta[,i] <- factor(meta[,i],levels = metaFactor[[i]])
    }
    ## mandatory annotation to be factor since in comparison ------
    if(!is.null(addFactor)){
        for(i in addFactor){
            if(!is.factor(meta[,i])){
                meta[,i] <- factor(meta[,i],levels=unique(meta[,i]))
            }
        }
    }
    return(meta)
}
metaFactor_checkNAempty <- function(meta){
    for(i in colnames(meta)){
        if(grepl("URL",i)) next
        if(is.character(meta[,i]) || is.factor(meta)){
            meta[,i] <- as.character(meta[,i])
            meta[is.na(meta[,i])|nchar(meta[,i])==0,i] <- "NA"
        }
    }
    return(meta)
}
metaFactor_addFactor <- function(meta,strMetaFactor,metaFactor=list()){
    for(i in colnames(meta)){
        if(grepl("URL",i)) next
        if(i %in%names(metaFactor)){
            if(sum(!meta[,i]%in%metaFactor[[i]])>0){
                metaFactor[[i]] <- c(metaFactor[[i]],unique(meta[!meta[,i]%in%metaFactor[[i]],i]))
            }
        }else{
            if(is.character(meta[,i]) || is.factor(meta))
                metaFactor[[i]] <- unique(meta[,i])
        }
    }
    metaFactor_saveYaml(metaFactor,strMetaFactor)
    return(metaFactor)
}
metaFactor_saveYaml <- function(ymlist,strF){
    cat(paste(paste0(names(ymlist),": ['",
                     sapply(ymlist,paste,collapse="','"),"']"),
              collapse="\n"),
        "\n",sep="",file=strF)
}
## EAqc functions ------
source("PC_Covariates.R")
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
plotPCanlaysis <- function(config,logTPM,meta,estC=NULL,effL=NULL){
    selCov <- unique(c(config$covariates_check,config$covariates_adjust))
    if(!is.null(selCov)) meta <- meta[,selCov]
    ## change the Well_Row from charactor to numeric
    oneMeta <- "Well_Row"
    if(oneMeta %in% colnames(meta)) meta[,oneMeta] <- as.numeric(as.factor(meta[,oneMeta]))
    # remove meta with only one unique value
    meta <- meta[,apply(meta,2,function(x)return(length(unique(x))>1)),drop=F]
    
    strPrefix <- paste0(config$output,"/notAdjusted")
    suppressMessages(suppressWarnings(
        Covariate_PC_Analysis(logTPM,meta,
                              out_prefix=strPrefix,
                              PC_cutoff=config$covariates_check_PCcutoff,
                              FDR_cutoff=config$covariates_check_FDRcutoff,
                              N_col=config$covariates_check_plotNcol)))
    message("-----> PC analysis without covariate adjusted:\n\t",strPrefix)
    
    if(!is.null(estC) && !is.null(effL)){
        if(is.null(config$covariates_adjust) || length(config$covariates_adjust)==0){
            warning("< covariates_adjust is NOT set in the config file, no covariate was adjusted! >")
        }else{
            message("====== removing covariates for visualization ...")
            batchX <- meta[,config$covariates_adjust,drop=F]
            logTPM <- suppressMessages(covariateRM(estC,effL,batchX=batchX,method='limma',
                                                   prior=config$count_prior))
            strPrefix <- paste0(config$output,"/Adjusted")
            suppressMessages(suppressWarnings(
                Covariate_PC_Analysis(logTPM,meta,
                                      out_prefix=strPrefix,
                                      PC_cutoff=config$covariates_check_PCcutoff,
                                      FDR_cutoff=config$covariates_check_FDRcutoff,
                                      N_col=config$covariates_check_plotNcol)))
            message("-----> PC analysis with covariate removal:\n\t",strPrefix)
            return(logTPM)
        }

    }
    return(NULL)
}
finishQC <- function(strMsg){
    message("==========================================")
    message("----->'EArun' can be used to obtain the QuickOmics object after necessary 'covariates_adjust' is set and comparison definition file is filled:")
    message("\t\t\t",strMsg$comparison_file)
    message("\t\tEArun ",strMsg$output,"/config.yml\n\n")
    message("-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.")
}

## EAsplit functions -------------
# EAdata: meta, counts, gInfo, effLength, logTPM, seqQC
splitUpdatePrj <- function(config,one,strOut){
    config$prj_name <- paste(config$prj_name,one,sep="_")
    config$prj_title <- paste0("(",one,") ",config$prj_name)
    config$output <- strOut
    config$DA_file_outpath <- paste0(strOut,"/DA_Import_Files")
    return(config)
}
splitSaveData <- function(X,strF,selRow=NULL,selCol=NULL,saveRowNames=NULL,...){
    message("\t",basename(strF))
    if(is.null(X)) return("")
    if(!is.null(selCol))
        X <- X[,selCol,drop=F]
    if(!is.null(selRow))
        X <- X[selRow,,drop=F]
    if(!is.null(saveRowNames))
        X <- setNames(cbind(rownames(X),X),c(saveRowNames,colnames(X)))
    write.table(X,strF,row.names=F,...)
    return(strF)
}
splitSaveFactor <- function(strSrc,strDest,rmOne){
    if(file.exists(strSrc)){
        cat(paste(grep(paste0("^",rmOne,":"),readLines(config$sample_factor),invert=T,value=T),
                  collapse="\n"),"\n",sep="",file=strDest)
    }
    return(strDest)
}
splitSaveComp <- function(strSrc,strDest,rmOne){
    comp <- read_file(strSrc,T)
    if(nrow(comp)>0){
        message("\textracting comparison information")
        selCom <- gsub(" ","",comp[,"Subsetting_group"])==paste(config$split_meta,rmOne,sep=":")
        if(sum(selCom)>0){
            comp <- comp[selCom,]
            comp[,"Subsetting_group"] <- ""
        }
    }
    return(splitSaveData(comp,strDest,saveRowNames="CompareName",sep=","))
}
splitSaveConfig <- function(config,strF){
    cat(paste(paste(names(config),
                sapply(config,function(x){
                    if(is.null(x)) return("")
                    if(length(x)==1) return(x)
                    return(paste0("['",paste(x,collapse="','"),"']"))
                }),
                sep=": "),collapse="\n"),
        "\n",sep="",file=strF)
    
}
splitOne <- function(config,EAdata,one){
    if(sum(EAdata$meta[,config$split_meta]==one)<2){
        message("ignore: ",one," with less than 2 samples")
        next
    }
    message("====== Creating sub project: ",one," ...")
    EAdata$meta <- EAdata$meta[EAdata$meta[,config$split_meta]==one,]
    strD <- paste0(config$output,"/",one,"/data")
    system(paste("mkdir -p",strD))
    strD <- normalizePath(strD)
    config <- splitUpdatePrj(config,one,dirname(strD))
    message("\t@",config$output)
    
    config$prj_counts <- splitSaveData(EAdata$counts,
                                       paste0(strD,"/count.tsv"),
                                       selCol=rownames(EAdata$meta),
                                       saveRowNames="gene_id",
                                       sep="\t")
    config$prj_counts <- splitSaveData(EAdata$effLength,
                                       paste0(strD,"/effLength.tsv"),
                                       selCol=rownames(EAdata$meta),
                                       saveRowNames="gene_id",
                                       sep="\t")
    config$prj_TPM <- splitSaveData(2^EAdata$logTPM-config$count_prior,
                                    paste0(strD,"/TPM.tsv"),
                                    selCol=rownames(EAdata$meta),
                                    saveRowNames="gene_id",
                                    sep="\t")
    config$prj_seqQC <- splitSaveData(EAdata$seqQC,
                                      paste0(strD,"/seqQC.tsv"),
                                      selRow=rownames(EAdata$meta),
                                      saveRowNames="Sample",
                                      sep="\t")
    config$gene_annotation <- splitSaveData(EAdata$gInfo,
                                      paste0(strD,"/geneAnnotation.csv"),
                                      saveRowNames="gene_id",
                                      sep=",")
    config$sample_meta <- splitSaveData(EAdata$meta,
                                        paste0(strD,"/sampleMeta.csv"),
                                        sep=",")
    config$sample_factor <- splitSaveFactor(config$sample_factor,
                                            paste0(strD,"/sampleMetaFactor.yml"),
                                            one)
    config$comparison_file <- splitSaveComp(config$comparison_file,
                                            paste0(strD,"/compareInfo.csv"),
                                            one)
    splitSaveConfig(config,paste0(dirname(strD),"/config.yml"))
    return(dirname(strD))
}
splitPrj <- function(config,EAdata){
    if(is.null(config$split_meta))
        stop("split_meta is not speicified in the config file")

    for(one in unique(EAdata$meta[,config$split_meta])){
        strOut <- splitOne(config,EAdata,one)
        suppressMessages(saveSessionInfo(paste0(strOut,"/session.EAsplit"),config$src))
    }
}
finishSplit <- function(){
    message("==========================================")
    message("Several ExpressionAnalysis project folders are created above")
    message("-----> 'EAqc' can be used to identify the covariates for each sub-project by using config in the folder:")
    message("-----> 'EArun' can be used to obtain the QuickOmics objects for each sub-project after comparison definition file is updated.\n")
}

## EArun functions ----
require(dplyr)
require(Hmisc)
source("QuickOmics_DEG.R")
saveCountsAlias <- function(config,estC){
    saveRDS(estC,file=paste0(config$output,"/",config$prj_name,"_estCount.rds"))
}
comparisonAnalysis <- function(config,estC,meta,comp_info){
    message("====== Starting DEG analyses ...")
    saveCountsAlias(config,estC)
    ## comparison -----------
    if(!is.null(config$qsub) && config$qsub){
        source(paste0(config$srcDir,"/qsubDEG.R"))
        return(qsubDEG(estC,meta,comp_info,config$output,config$srcDir,core=config$core))
    }else{
        return(Batch_DEG(estC,meta,comp_info,core=config$core))
    }
}
Hmisc.rcorr <- function (x, y, type = "pearson"){
    type <- match.arg(type)#c("pearson", "spearman")
    if (!missing(y))
        x <- cbind(x, y)
    x[is.na(x)] <- 1e+50
    storage.mode(x) <- "double"
    p <- as.integer(ncol(x))
    if (p < 1)
        stop("must have >1 column")
    n <- as.integer(nrow(x))
    if (n < 5)
        stop("must have >4 observations")
    h <- .Fortran(Hmisc:::F_rcorr, x, n, p, itype = as.integer(1 + (type =="spearman")),
                  hmatrix = double(p * p), npair = integer(p * p), double(n),
                  double(n), double(n), double(n), double(n),
                  integer(n))
    npair <- matrix(h$npair, ncol = p)
    h <- matrix(h$hmatrix, ncol = p)
    h[h > 1e+49] <- NA
    nam <- dimnames(x)[[2]]
    dimnames(h) <- list(nam, nam)
    dimnames(npair) <- list(nam, nam)
    #https://stackoverflow.com/questions/63994852/r-in-sqrt1-h-h-nans-produced-from-within-rcorr-full-sample-data-avai
    #P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/sqrt(1 - h * h), npair - 2)), ncol = p)
    P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/max(0, 1-h^2), npair - 2)), ncol = p)
    P[abs(h) == 1] <- 0
    diag(P) <- NA
    dimnames(P) <- list(nam, nam)
    structure(list(r = h, n = npair, P = P), class = "rcorr")
}
saveNetwork <- function(X,config){
    message("Obtaining networks ...")
    cor_cutoff <- config$gene_network_cor_cutoff
    p_cutoff <- config$gene_network_p_cutoff
    variableN <- config$gene_network_high_variable_N
    edge_max <- as.numeric(config$gene_network_max_edge)
    edge_min <- as.numeric(config$gene_network_min_edge)

    suppressMessages(require(tibble))
    if(nrow(X)>variableN){
        X <- X[order(apply(X,1,sd),decreasing=T)[1:variableN],]
    }
    print(system.time(cor_res <- Hmisc.rcorr(as.matrix(t(X)))))
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
    network <- network %>% dplyr::mutate_if(is.factor, as.character) %>%
        dplyr::filter(!is.na(cor) & abs(cor) > cor_cutoff & p < p_cutoff)
    save(network,file=paste0(config$output,"/",config$prj_name,"_network.RData"))
}
updateMeta <- function(config,meta){
    comp <- read_file(config$comparison_file,T)
    meta <- metaFactor(meta,config$sample_factor,unique(comp$Group_name))
    if(!"group" %in% colnames(meta)){
        selG <- c(grep("^group$",colnames(meta),ignore.case=T,value=T),
                  comp[1,"Group_name"])[1]
        meta <- cbind(meta,group=meta[,selG])
    }
    return(meta)
}
formatQuickOmicsResult <- function(DEGs,logTPM,grp,gInfo){
    Dw <- data.frame(gInfo[rownames(logTPM),c("UniqueID","Gene.Name",'id')],
                     Intensity=apply(logTPM,1,mean))
    for(i in unique(grp)){
        if(sum(grp==i)<1){
            next
            message("===== warning: no sample for ",i)
        }else if(sum(grp==i)==1){
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
                      meta[,-grep("^group$",colnames(meta)),drop=F])
    suppressWarnings(MetaData[is.na(MetaData)] <- "")
    return(MetaData)
}
saveQuickOmics <- function(config,EAdata,DEGs){
    message("saving QuickOmics object ...")
    EAdata$meta <- updateMeta(config,EAdata$meta)
    comp_info <- read_file(config$comparison_file,T)
    
    message("\tFormating the expression data")
    data_wide <- EAdata$logTPM
    data_long <- melt(as.matrix(EAdata$logTPM))
    colnames(data_long) <- c("UniqueID","sampleid","expr")
    config$sample_group <- ifelse(length(config$sample_group)<1,"group",config$sample_group)
    data_long <- cbind(data_long,group=EAdata$meta[data_long$sampleid,config$sample_group])
    
    message("\tFormating the DEG results")
    compRes <- formatQuickOmicsResult(DEGs,EAdata$logTPM,EAdata$meta[,"group"],EAdata$gInfo)
    data_results <- compRes$Dw
    results_long <- compRes$Dl
    
    message("\tFormating the sample meta information")
    MetaData <- formatQuickOmicsMeta(EAdata$meta,names(DEGs))
    ProteinGeneName <- EAdata$gInfo
    
    message("\tsaving...")
    save(data_results,results_long,
         data_wide,data_long,
         MetaData,ProteinGeneName,
         comp_info,
         file=paste0(config$output,"/",config$prj_name,".RData"))
    ## save the project csv file -------
    write.csv(data.frame(Name=ifelse(is.null(config$prj_title),
                                     config$prj_name,
                                     paste(config$prj_name,config$prj_title,sep=": ")),
                         ShortName=config$prj_name,
                         ProjectID=config$prj_name,
                         Species=config$species, 
                         ExpressionUnit=config$ylab,
                         Path=config$output),
              file=paste0(config$output,"/",config$prj_name,".csv"),
              row.names=F)
}
finishRun <- function(strMsg){
    message("=================================================\nResults are saved in ",strMsg$output)
    system(paste0("cp ",strMsg$output,"/",strMsg$prj_name,"* ",strMsg$QuickOmics_test_folder))
    message(paste0("\n-----> Please visit: ",strMsg$QuickOmics_test_link,strMsg$prj_name))
    message("Please carefully review the results before publishing:")
    message("----->'EApub' can be used to publish the project into ShinyOne project manager: ",
            sapply(strsplit(strMsg$shinyApp,"\\/"),
                   function(x)return(paste(x[1:grep("shiny",x)],collapse="/"))))
    
}


## EApub functions --------
getShinyOneInfo <- function(config){
    shinyOneData <- config[grep("^shinyOne_",names(config))]
    names(shinyOneData) <- gsub("shinyOne_","",names(shinyOneData))
    
    shinyOneData[['Title']] <- ifelse(is.null(config$shinyOne_Title),
                                      config$prj_name,
                                      config$shinyOne_Title)
    shinyOneData[["Species"]] <- config$species
    shinyOneData[["TSTID"]] <- ifelse(grepl("^TST",config$prj_name),substr(config$prj_name,1,8),"")
    
    shinyOneData[["Data_Generated_By"]] <- paste(system("whoami",intern=T),
                                                 shinyOneData[["Data_Generated_By"]],
                                                 sep="; ")
    shinyOneData[["Date"]] <- as.character(Sys.Date())
    shinyOneData[["URL"]] <- paste0(config$QuiclOmics_publish_link,config$prj_name)
    return(shinyOneData)
}
pubShinyOne <- function(config){
    strF <- paste0(config$QuickOmics_publish_folder,config$prj_name,".RData")
    if(file.exists(strF)){
        stop("The project already exists in ShinyOne!\nPlease remove the record and associated files or change prj_name and re-EArun!")
    }
    message("preparing information for ShinyOne")
    shinyOneData <- getShinyOneInfo(config)
    shinyOneCMD <- paste0("curl -s -k -X POST -d 'data={",
                          paste(paste0('"',names(shinyOneData),'": "',shinyOneData,'"'),collapse = ", "),
                          "}' '",config$shinyApp,"api_add_project.php?api_key=lnpJMJ5ClbuHCylWqfBY8BoxxdrpU0'")

    message("submitting to ShinyOne manager ...")
    res <- system(shinyOneCMD,intern=T)
    shinyMsg <- tryCatch({
        rjson::fromJSON(res)
    },error=function(eMsg){
        stop(paste0(paste(res,collapse="\n"),
                    "\nPlease contact Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]"))
    })
    if(!shinyMsg$Status){
        stop(paste0(paste(paste(names(shinyMsg),shinyMsg,sep=":"),collapse="\n"),
                    "\nPlease contact Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]"))
    }
    system(paste0("cp ",config$output,"/",config$prj_name,"* ",config$QuickOmics_publish_folder))
    return(shinyMsg$ID)
}
finishShinyOne <- function(shinyMsg){
    message("=================================================\nShinyOne access: ",
            shinyMsg$shinyApp,"app_project_review.php?ID=",shinyMsg$ID)
    
}

## others ----------
