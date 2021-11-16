#source("utility.R",chdir=T)

.libPaths(c(grep("home",.libPaths(),invert=T,value=T),grep("home",.libPaths(),value=T)))
require(data.table)

## functions used for EAinit -----
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
        pInfo[["gAnnotion"]] <- getAnnotation(strConfig,strGannotation)
        
        strPrj <- gsub("tsv$","json",strSample)
        if(file.exists(strPrj)){
            prjInfo <- rjson::fromJSON(file=strPrj)$Project
            pInfo[["prjID"]] <- prjInfo$TSTID
            pInfo[["prjTitle"]] <- prjInfo$Study_Title
        }
        
        pInfo[["strCount"]] <- normalizePath(paste0(strInput,"/combine_rsem_outputs/genes.estcount_table.txt"))
        pInfo[["strEffLength"]] <- normalizePath(getEffectLength(strInput))
        pInfo[["strSeqQC"]] <- normalizePath(paste0(strInput,"/combine_rnaseqc/combined.metrics.tsv"))
    }
    return(pInfo)
}
getEffectLength<-function(strPath){
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
getAnnotation <- function(strF,strPath=NULL){
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
cleanTST <- function(strSrc,strDest=NULL,sep="\t",beReturn=F){
    D <- read.table(strSrc,sep="\t",header=T,as.is=T,check.names=F,quote="")
    if(sum(grepl("^TST",colnames(D)))>0){
        colnames(D) <- sapply(strsplit(sapply(strsplit(colnames(D),"\\|"),head,1),
                                       "_"),
                              function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
    }else if(sum(grepl("^TST",D[,1]))>0){
        D[,1] <- sapply(strsplit(gsub(".genome.sorted","",D[,1]),
                                 "_"),
                        function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
    }
    if(!is.null(strDest)) write.table(D,strDest,row.names=F,sep="\t")
    if(beReturn) return(D)
}
appendMeta <- function(pInfo,sample_name,selQC){
    message("Appending sequencing QC into sample meta file ...")
    if(!is.null(pInfo$sInfo) && !is.null(pInfo$strSeqQC)){
        qc <- cleanTST(pInfo$strSeqQC,beReturn=T)
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
        
        configTmp <- gsub("^shinyOne_Title:",paste0("shinyOne_Title: ",pInfo$prjID,"-",pInfo$prjTitle),configTmp)
        
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

## general functions used by all ------
checkConfig <- function(config){
    if(is.null(config$sample_name))
        stop("sample_name is required. Default is Sample_Name")
    
}
saveSeesionInfo <- function(strF,strSRC){
    sink(strF)
    cat("EA version: git-",
        system(paste('git --git-dir',normalizePath(paste0(strSRC,"../.git")),'rev-parse HEAD'),
               intern=T),sep="","\n\n")
    print(sessionInfo())
    sink()
    
}



