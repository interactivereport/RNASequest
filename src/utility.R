strConfig <- file.path(getwd(),"sys.yml")
if(!file.exists(strConfig))
  stop(paste0("=====\nPlease contact admins to set the sys.yml in ",args[1],".\nAn Example is 'sys_example.yml'.\n====="))
## save session ------
initialMsg <- function(strSRC="."){
  message("\n\n***** ",Sys.time()," *****")
  strSRC <- normalizePath(strSRC)
  if(dir.exists(file.path(strSRC,".git"))){
    gitConfig <- readLines(file.path(strSRC,".git","config"))
    pos <- grep("^\\[.*\\]$",gitConfig)
    sel <- (1:length(pos))[grepl("remote",gitConfig[pos])&grepl("origin",gitConfig[pos])]
    if(length(sel)==0) return()
    pos <- c(pos,length(gitConfig)+1)
    url <- sapply(strsplit(grep("^url",trimws(gitConfig[pos[sel]:(pos[sel+1]-1)]),value=T),
                           " "),tail,1)
    gitLog <- unlist(strsplit(unlist(tail(data.table::fread(file.path(strSRC,".git","logs","HEAD"),header=F),1))[1]," "))
    message("###########\n## ExpressionAnalysis: ",url)
    message("## Pipeline Path: ",strSRC)
    message("## Pipeline Date: ",
            format(as.POSIXct(as.numeric(tail(gitLog,2)[1]),
                              origin="1970-01-01"),
                   format="%Y-%m-%d %H:%M:%S"),
            " ",tail(gitLog,1))
    message("## git HEAD: ",gitLog[2],"\n###########\n")
  }else if(file.exists(file.path(strSRC,"release"))){
    system(paste("cat",file.path(strSRC,"release")))
    message("## Pipeline Path: ",strSRC)
  }
}
saveSessionInfo <- function(strF,strSRC){
    message("\nPowered by ",sysConfig$powerby)
    writeLines(capture.output(sessionInfo()), strF)
    conn <- file(strF,"a")
    sink(conn,type="message")
    initialMsg(strSRC)
    sink(type="message")
    close(conn)
}
## EAinit functions -----
source("gtf.gz_to_gene_info.R")
checkInputDir <- function(strInput,sysConfig=NULL){
    strInput <- normalizePath(strInput)
    if(!dir.exists(strInput)){
        stop(paste(strInput, "is not a valid path"))
    }
    #for internal data structure
    prjInfo <- initInternal(strInput,sysConfig)
    if(is.null(prjInfo)){
        prjInfo <- initFileNameMatch(strInput,sysConfig)
    }
    return(prjInfo)
}

initInternal <- function(strInput,sysConfig){
    pInfo <- NULL
    strSample <- paste0(strInput,"/samplesheet.tsv")
    strConfig <- list.files(strInput,"^analysis-",full.names=T)
    if(file.exists(strSample) && length(strConfig)>0){
        pInfo <- list()
        pInfo[["sInfo"]] <- data.frame(fread(strSample, header=T),stringsAsFactors=F)
        pInfo$sInfo <- pInfo$sInfo[,apply(pInfo$sInfo,2,function(x)return(sum(!is.na(x))>0))]
        strConfig <- strConfig[1]
        message("\tusing ",strConfig)
        analysisP <- yaml::read_yaml(strConfig)
        if('properties' %in% names(analysisP)){
            analysisP <- analysisP$properties
        }
        pInfo[["prjID"]] <- pInfo[["prjTitle"]] <- getAnalysisInfo(analysisP,sysConfig$DNAnexus$prjID) #analysisP[[sysConfig$DNAnexus$prjID]]
        pInfo[["uName"]] <- getAnalysisInfo(analysisP,sysConfig$DNAnexus$uName) #analysisP[[sysConfig$DNAnexus$uName]]
        pInfo[["species"]] <- getAnalysisInfo(analysisP,sysConfig$DNAnexus$species) #analysisP[[sysConfig$DNAnexus$species]]
        pInfo[["gAnnotion"]] <- extractAnnotation(getAnalysisInfo(analysisP,sysConfig$DNAnexus$species), #analysisP[[sysConfig$DNAnexus$species]],
        										  getAnalysisInfo(analysisP,sysConfig$DNAnexus$ref),#analysisP[[sysConfig$DNAnexus$ref]],
                                                  sysConfig$genome_path)
        
        strPrj <- gsub("tsv$","json",strSample)
        if(file.exists(strPrj)){
            prjInfo <- rjson::fromJSON(file=strPrj)$Project
            if(!is.null(prjInfo$Study_Title))
              pInfo[["prjTitle"]] <- prjInfo$Study_Title
        }
        #res <- list.files(strInput,"estcount",recursive=T,full.names=T)
        
        
        pInfo[["strCount"]] <- getDNAnexusFile(strInput,sysConfig$DNAnexus$count,"count")
        pInfo[["strEffLength"]] <- getEffectLengthFile(strInput,
                                                       sysConfig$DNAnexus$effL,
                                                       sysConfig$DNAnexus$indFlag)
        if(is.null(pInfo[["strEffLength"]])){
            pInfo[["strTPM"]] <- getDNAnexusFile(strInput,sysConfig$DNAnexus$tpm,"TPM")
        }
        pInfo[["strSeqQC"]] <- getDNAnexusFile(strInput,sysConfig$DNAnexus$seqQC,"seqQC")
        pInfo <- getCovariates(pInfo,sysConfig$notCovariates)
        pInfo$datatype <- "DNAnexus"
    }
    return(pInfo)
}
getAnalysisInfo <- function(analysisInfo,keys){
	aInfo <- NULL
	for(one in keys){
		aInfo <- analysisInfo[[one]]
		if(!is.null(aInfo)) return(aInfo)
	}
	stop(paste("None of the pattern can be found in analysis information:",paste(one,collapse=";")))
}
extractAnnotation <- function(species,sVersion,genome_path=NULL){
    message("Create gene annotation ...")
    gtfPath <- list.files(paste0(genome_path,"/rnaseq/",
                                 species,"/",sVersion),
                          "gtf.gz$",full.names=T)[1]
    if(is.na(gtfPath)) stop(paste("missing reference GTF for",species,sVersion))
    strF <- gsub("gtf.gz$","gene_info.csv",gtfPath)
    if(!file.exists(strF)){
    	message("\t",gtfPath)
        gtf.gz_to_gene_info(gtfPath)
    }

    #columns in strF: unique ID, gene_name, gene_type,....
    message("\t",strF)
    gInfo <- read.csv(strF,as.is=T,check.names=F)
    gInfoName <- colnames(gInfo)
    gInfo <- cbind(0:(nrow(gInfo)-1),gInfo)
    dimnames(gInfo) <- list(gInfo[,2],
                            c("id","UniqueID","Gene.Name","Biotype",tail(gInfoName,-3)))
    ## if the lookup table is available
    strLookup <- paste0(genome_path,"/rnaseq/lookup/",
                        sVersion,".csv")
    if(file.exists(strLookup)){
        gLookup <- read.csv(strLookup,row.names=1,as.is=T,check.names = F)
        gInfo <- merge(gInfo,gLookup,by="row.names",all=T,sort=F)
        rownames(gInfo) <- gInfo[,1]
        gInfo <- gInfo[,-1]
    }
    return(gInfo)
}
getDNAnexusFile <- function(strInput,pattern,key){
    message("getting ",key," file by ",paste(pattern,collapse=" or "))
	res <- NULL
	for(one in pattern){
		res <- list.files(strInput,one,full.names=T,recursive =T)
		if(length(res)>0){
			message("\t",res[1])
			return(normalizePath(res[1]))
		}
	}
	stop("No ",key," matrix found")
    #res <- list.files(strInput,pattern,full.names=T)
    #if(length(res)==0){
    #    res <- list.files(paste0(strInput,"/combine_rsem_outputs"),
    #                      pattern,full.names=T)
    #}
    
    #return(normalizePath(res[order(nchar(res))][1]))
}
getEffectLengthFile <-function(strInput,pattern,indFlag){
    message("Extracting effective length by ",paste(pattern,collapse=" or "))
	for(one in pattern){
		res <- list.files(strInput,pattern,full.names=T,recursive =T)
		if(length(res)>0){
			message("\t",res[1])
			return(normalizePath(res[1]))
		}
	}
    ## extract from each individual files
    strFs <- list.files(strInput,indFlag,full.names=T)
    if(length(strFs)==0){
        strFs <- list.files(paste0(strInput,"/rsem"),indFlag,full.names=T)
    }
    if(length(strFs)==0){
        message("No effective length, TPM is used for visualization (covariate adjustment is disabled).")
        return(NULL)
    }
    
    D <- NULL
    selColumn <- c("gene_id","effective_length")
    for(i in strFs){
        X <- fread(i,sep="\t",header=T)[,..selColumn]#,as.is=T,row.names=1
        colnames(X) <- c(selColumn[1],gsub(paste0(indFlag,"|\\.",indFlag),paste0("|",selColumn[2]),basename(i)))
        if(is.null(D)) D <- X
        else D <- merge(D,X,by=selColumn[1],all=T)
    }
    strF <- paste0(strInput,"/",pattern,".tsv")
    write.table(D,file=strF,sep="\t",row.names=F,quote=F)
    return(strF)
}
getQCfile <- function(strInput,pattern){
    message("\tgetting ",pattern)
	for(one in pattern){
		res <- list.files(strInput,pattern,full.names=T,recursive =T)
		if(length(res)>0){
			message("\t\t",res[1])
			return(normalizePath(res[1]))
		}
	}
    res <- list.files(strInput,pattern,full.names=T)
    if(length(res)==0){
        res <- list.files(paste0(strInput,"/combine_rnaseqc"),
                          pattern,full.names=T)
    }
    if(length(res)==0) stop("Internal alignment QC file is missing!")
    return(normalizePath(res[order(nchar(res))][1]))
}

initFileNameMatch <- function(strInput,sysConfig){
    pInfo <- NULL
    fList <- list.files(strInput)
    strF <- setNames(rep("",length(sysConfig$FileName)),names(sysConfig$FileName))
    for(one in names(strF)){
        if(length(sysConfig$FileName[[one]])>1) next
        a <- grep(sysConfig$FileName[[one]],fList,value=T,ignore.case = T)
        strF[one] <- ifelse(length(a)==0,"",normalizePath(paste0(strInput,"/",a[1])))
    }
    if(file.exists(strF["prj_counts"]) && file.exists(strF["sample_meta"]) &&
       (file.exists(strF["prj_TPM"]) || file.exists(strF["prj_effLength"]))){
        pInfo <- list(prjID=basename(strInput),
                      prjTitle=basename(strInput),
                      species="human")
        if(file.exists(strF["projectFile"])){
            prjInfo <- yaml::read_yaml(strF["projectFile"])
            for(one in names(sysConfig$FileName$projectEntry)){
                prjInfo[[one]] <- prjInfo[[sysConfig$FileName$projectEntry[[one]]]]
            }
        }
        if(file.exists(strF["gene_annotation"])){
            pInfo[["gAnnotion"]] <- data.frame(fread(strF["gene_annotation"], header=T),stringsAsFactors=F)
            rownames(pInfo[["gAnnotion"]]) <- pInfo[["gAnnotion"]][,1]
            pInfo[["gAnnotion"]] <- pInfo[["gAnnotion"]][,-1]
        }
        pInfo[["strCount"]] <- strF["prj_counts"]
        
        if(file.exists(strF["prj_effLength"]))
            pInfo[["strEffLength"]] <- strF["prj_effLength"]
        if(file.exists(strF["prj_TPM"]))
            pInfo[["strTPM"]] <- strF["prj_TPM"]
        pInfo[["sInfo"]] <- data.frame(fread(strF["sample_meta"], header=T),stringsAsFactors=F)
        colnames(pInfo[["sInfo"]])[1] <- "Sample_Name"
        if(file.exists(strF["prj_seqQC"]))
            pInfo[["strSeqQC"]] <- strF["prj_seqQC"]
        pInfo <- getCovariates(pInfo,sysConfig$notCovariates)
        if(file.exists(strF[["comparison_file"]]))
            pInfo[["comparison_file"]] <- data.frame(fread(strF[["comparison_file"]], header=T),stringsAsFactors=F)
        
        pInfo$datatype <- "FileMatch"
    }
    return(pInfo)
}

createComparison <- function(strComp,comp=NULL){
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
    if(is.null(comp)){
        message("Create empty comparison template ...")
        cat(paste(comTitle,collapse=","),"\n",sep="",file=strComp)
    }else{
        compTable <- matrix("",nrow=nrow(comp),ncol=length(comTitle),dimnames=list(1:nrow(comp),comTitle))
        compTable[,"CompareName"] <- apply(comp,1,function(x)return(paste0(x[1],"_",x[2],".vs.",x[3])))
        compTable[,"Group_name"] <- comp[,1]
        compTable[,"Group_test"] <- comp[,2]
        compTable[,"Group_ctrl"] <- comp[,3]
        compTable[,"Analysis_method"] <- "DESeq2"
        write.csv(compTable,file=strComp,row.names=F)
    }
}
saveGeneAnnotation <- function(gAnno,strF){
    if(!is.null(gAnno)){
        write.csv(gAnno,file=strF)
    }
}
cleanTST <- function(strSrc,strDest=NULL,sep="\t"){
    if(is.null(strSrc) || !file.exists(strSrc)) return(NULL)
    D <- read.table(strSrc,sep=sep,header=T,as.is=T,check.names=F,quote="")
    if(nchar(colnames(D)[1])==0) colnames(D)[1] <- "gene_id"
    if(sum(grepl("^TST",colnames(D)))>0){
        colnames(D) <- sapply(strsplit(sapply(strsplit(colnames(D),"\\|"),head,1),
                                       "_"),
                              function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
    }else if(sum(grepl("^TST.*.genome.sorted$",D[,1]))>0){
        D[,1] <- sapply(strsplit(gsub(".genome.sorted","",D[,1]),
                                 "_"),
                        function(x)return(paste(grep("^TST",x,invert=T,value=T),collapse="_")))
        D <- D[!duplicated(D[,1])&nchar(D[,1])>0,]
    }else{
    	#fastr 3
    	res <- do.call(rbind.data.frame,strsplit(colnames(D),"\\|"))
    	sel <- apply(res,2,function(x)return(length(unique(x[-1]))))>1
    	if(sum(!sel)>0) colnames(D) <- apply(res[,sel,drop=F],1,paste,collapse="|")
    }
    if(!is.null(strDest)) write.table(D,strDest,row.names=F,sep=sep)
    else return(D)
}
appendMeta <- function(pInfo,sample_name,selQC){
    message("Appending sequencing QC into sample meta file ...")
    if(!is.null(pInfo$sInfo) && !is.null(pInfo$strSeqQC)){
        qc <- cleanTST(pInfo$strSeqQC)
        rownames(qc) <- qc[,1]
        if(sum(!selQC%in%colnames(qc))>0){
          message("\tWarning: Missing sequence QC: ",paste(selQC[!selQC%in%colnames(qc)],collapse=", "))
          selQC <- selQC[selQC%in%colnames(qc)]
        }
        #selQC <- selQC[selQC%in%colnames(qc)]
        qc <- qc[,selQC,drop=F]
        colnames(qc) <- gsub("^3","x3",gsub("5","x5",gsub("'","p",gsub(" ","_",gsub("\\%","percentage",colnames(qc))))))
        ## only the sequenced samples will be included (remove not sequenced samples from sample sheet)
        message(sample_name)
        print(head(pInfo$sInfo))
        pInfo$sInfo <- pInfo$sInfo[pInfo$sInfo[,sample_name]%in%rownames(qc),]
        
        
        rownames(pInfo$sInfo) <- pInfo$sInfo[,sample_name]
        
        #if(sum(!rownames(pInfo$sInfo)%in%rownames(qc))>0)
        #    stop(paste0("Sample Names (",
        #               paste(rownames(pInfo$sInfo)[!rownames(pInfo$sInfo)%in%rownames(qc)],collapse=", "),
        #               ") specified in sample sheet cannot be found in seqQC file ",
        #               pInfo$strSeqQC))
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
createInit <- function(strInput,configTmp,pInfo,sysConfig=NULL){
    strInput <- normalizePath(strInput)
    message("Creating project folder")
    ix <- 0
    while(dir.exists(strOut<-paste0(strInput,"/EA",gsub("\\-","",Sys.Date()),"_",ix))){
        ix <- ix+1
    }
    system(paste("mkdir -p",paste0(strOut,"/data")))
    
    strCount <- paste0(strOut,"/data/count.tsv")
    strEffLength <- paste0(strOut,"/data/effLength.tsv")
    strTPM <- paste0(strOut,"/data/TPM.tsv")
    strSeqQC <- paste0(strOut,"/data/seqQC.tsv")
    
    strMeta <- paste0(strOut,"/data/sampleMeta.csv")
    strMetaFactor <- paste0(strOut,"/data/sampleMetaFactor.yml")
    strGinfo <- paste0(strOut,"/data/geneAnnotation.csv")
    saveGeneAnnotation(pInfo$gAnnotion,strGinfo)
    
    strComp <- paste0(strOut,"/data/compareInfo.csv")
    createComparison(strComp,pInfo[["comparison_file"]])

    message("saving initialization ...")
    if(is.null(pInfo)){
        configTmp <- gsub("initPrjName"," #required",configTmp)
        configTmp <- gsub("initPrjTitle"," #optional",configTmp)
        configTmp <- gsub("initCounts"," #required",configTmp)
        configTmp <- gsub("initEffLength"," #if provided prj_TPM is ignored",configTmp)
        configTmp <- gsub("initSeqQC"," #optional",configTmp)
        configTmp <- gsub("initTPM"," #prj_effLength or prj_TPM is required",configTmp)
        configTmp <- gsub("initPrjMeta"," #required",configTmp)
        configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
        configTmp <- gsub("initSpecies"," #required",configTmp)
        configTmp <- gsub("initGeneAnnotation"," #optional",configTmp)
        configTmp <- gsub("initOutput",strOut,configTmp)
        configTmp <- gsub("initCovariates","",configTmp)
        configTmp <- gsub("initPrjComp",strComp,configTmp)
    }else if(pInfo$datatype=="DNAnexus"){
        cleanTST(pInfo[["strCount"]],strDest=strCount)
        if(is.null(pInfo[["strEffLength"]])){
            cleanTST(pInfo[["strTPM"]],strDest=strTPM)
            configTmp <- gsub("initTPM",strTPM,configTmp)
            configTmp <- gsub("initEffLength","",configTmp)
        }else{
            cleanTST(pInfo[["strEffLength"]],strDest=strEffLength)
            configTmp <- gsub("initTPM","",configTmp)
            configTmp <- gsub("initEffLength",strEffLength,configTmp)
        }
        cleanTST(pInfo[["strSeqQC"]],strDest=strSeqQC)
        write.csv(pInfo$sInfo,file=strMeta,row.names=F)
        
        configTmp <- gsub("initPrjName",pInfo[["prjID"]],configTmp)
        configTmp <- gsub("initPrjTitle",pInfo[["prjTitle"]],configTmp)
        configTmp <- gsub("initCounts",strCount,configTmp)
        configTmp <- gsub("initEffLength",strEffLength,configTmp)
        configTmp <- gsub("initSeqQC",strSeqQC,configTmp)
        
        
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
                          paste("shinyOne_Data_Generated_By:",pInfo$uName),
                          configTmp)
    }else if(pInfo$datatype=="FileMatch"){
        file.copy(pInfo[["strCount"]],strCount)
        write.csv(pInfo$sInfo,file=strMeta,row.names=F)
        
        configTmp <- gsub("initPrjName",pInfo[["prjID"]],configTmp)
        configTmp <- gsub("initPrjTitle",pInfo[["prjTitle"]],configTmp)
        configTmp <- configTmp <- gsub("initCounts",strCount,configTmp)
        if(!is.null(pInfo[["strEffLength"]]) && file.exists(pInfo[["strEffLength"]])){
            file.copy(pInfo[["strEffLength"]],strEffLength)
            configTmp <- gsub("initEffLength",strEffLength,configTmp)
        }else{
            configTmp <- gsub("initEffLength"," #if provided prj_TPM is ignored",configTmp)
        }
        if(!is.null(pInfo[["strTPM"]]) && file.exists(pInfo[["strTPM"]])){
            file.copy(pInfo[["strTPM"]],strTPM)
            configTmp <- gsub("initTPM",strTPM,configTmp)
        }else{
            configTmp <- gsub("initTPM"," #prj_effLength or prj_TPM is required",configTmp)
        }
        if(!is.null(pInfo[["strSeqQC"]]) && file.exists(pInfo[["strSeqQC"]])){
            file.copy(pInfo[["strSeqQC"]],strSeqQC)
            configTmp <- gsub("initSeqQC",strSeqQC,configTmp)
        }else{
            configTmp <- gsub("initSeqQC"," #optional",configTmp)
        }
        configTmp <- gsub("initPrjMeta",strMeta,configTmp)
        configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
        configTmp <- gsub("initSpecies",pInfo$species,configTmp)
        if(!is.null(pInfo$gAnnotion)){
            configTmp <- gsub("initGeneAnnotation",strGinfo,configTmp)
        }else{
            configTmp <- gsub("initGeneAnnotation"," #optional",configTmp)
        }
        configTmp <- gsub("initOutput",strOut,configTmp)
        configTmp <- gsub("initCovariates",paste0("[",paste(pInfo$covariates,collapse=","),"]"),configTmp)
        configTmp <- gsub("initPrjComp",strComp,configTmp)
        
    }else{
        stop("Unknown data input")
    }
    ## cellmap availability
    if(!is.null(cm_profiles <- cellmap_avaid_profiles(sysConfig))){
        configTmp <- c(configTmp,"# ====== This section for CellMap is available to estimate cell type proportions =========")
        configTmp <- c(configTmp,paste0("CellMap_profile:  #please specify one profiler (human) from ",paste(cm_profiles,collapse=", ")))
        configTmp <- c(configTmp,paste0("CellMap_rm_ct:  #please specify the cell type(s) to be removed from above profile"))
        configTmp <- c(configTmp,paste0("CellMap_species:  #please specify the species, human/mouse supported"))
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
    correctNaive(config)
    if(is.null(config$prj_name) || nchar(config$prj_name)<2){
      stop("'prj_name' is required in the config.")
    }
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
        if(is.null(config$covariates_method))
            config$covariates_method <- "limma"
        if(!config$covariates_method%in%c("limma","combat"))
            stop(paste("Unknown covariate adjust method:",config$covariates_method))
        if(grepl("combat",config$covariates_method) && length(config$covariates_adjust)>1)
            stop("Combat only supports 1 covarite adjustment!")
    }
    return(config)
}
getEAData <- function(config,withCom=F){
    D <- getMeta(config)
    if(withCom)
        D$comp_info <- checkComparisonInfo(read_file(config$comparison_file,T),
                                           D$meta,config$comparison_file)
    D <- c(D,getCounts(config,rownames(D$meta)))
    D <- c(D,getEffLength(config,colnames(D$counts),rownames(D$counts),D$gInfo,config$min_median_effective_length))
    D <- c(D,getTPM(config,colnames(D$counts),rownames(D$counts)))
    D <- c(D,getSeqQC(config,colnames(D$counts)))
    if(is.null(D$logTPM)){
        D$logTPM <- as.data.frame(covariateRM(D$counts,D$effLength,prior=config$count_prior))
    }
    if(is.null(D$gInfo)){
        D$gInfo <- data.frame(row.names = rownames(D$counts),
                              id=0:(nrow(D$counts)-1),
                              UniqueID=rownames(D$counts),
                              Gene.Name=rownames(D$counts))
    }
    return(filterGene(config,D))
}
getMeta <- function(config){
    message("reading sample meta")
    meta <- as.data.frame(data.table::fread(config$sample_meta))
    checkMeta(meta,config)
    rownames(meta) <- meta[,config$sample_name]
    return(setMetaFactor(meta,config$sample_factor))
}
checkMeta <- function(meta,config){
    message("\t",nrow(meta)," samples in sample meta file")
    stopifnot(nrow(meta)>1)
    message("\tchecking against config file")
    if(!config$sample_name%in%colnames(meta))
        stop(paste0("sample_name (",config$sample_name,
                    ") is NOT a column in sample meta file (",
                    config$sample_meta,")!"))
    if(sum(duplicated(meta[,config$sample_name])))
        stop(paste0("sample_name column (",config$sample_name,") contains duplicates in sample meta file."))
    if(!is.null(config$sample_alias)){
        if(!config$sample_alias%in%colnames(meta))
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
    #D <- read.table(config$prj_counts,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    D <- data.table::fread(config$prj_counts,header=T,check.names=F)
    D <- data.frame(row.names=D[[1]],D[,-1],check.names=F)
    
    if(ncol(D)<2) stop("No sample detected!")
    ix <- apply(as.matrix(D),1,function(x)return(sum(x>=config$min_count)))>=config$min_sample
    message("\tFiltering genes (",sum(ix),") with minimal counts (>=) ",
            config$min_count," in at least (>=) ",config$min_sample," samples")
    D <- D[ix,]
    gInfo <- NULL
    if(!is.null(config$gene_annotation)){
        gInfo <-  read.csv(config$gene_annotation,row.names=1,as.is=T)
        if(sum(!c("UniqueID","Gene.Name",'id')%in%colnames(gInfo))>0)
            stop("Columns of 'UniqueID', 'Gene.Name' and 'id' have to be defined in gene annotation file")
        gID <- intersect(rownames(D),rownames(gInfo))
        if(length(gID)<nrow(D)){
            message("\tFiltering genes (",length(gID),") from (",nrow(D),") by gene annotation")
            D <- D[gID,]
            gInfo <- gInfo[gID,]
        }
    }
    D <- checkSampleName(D,sID)
    return(filterGene(config,list(counts=D,gInfo=gInfo)))
}
checkSampleName <- function(D,sID){
    if(sum(!sID%in%colnames(D))>0)
        stop(paste0("samples (",paste(c(head(sID[!sID%in%colnames(D)],5),"..."),collapse=","),
                    ") defined in sample meta table are NOT available in count matrix"))
    delSample <- colnames(D)[!colnames(D)%in%sID]
    if(length(delSample)>0){
        message("The following ",length(delSample)," samples are removed from count table since they are not defined in sample meta:")
        message("\t",paste(head(delSample,5),collapse="; "),ifelse(length(delSample)>5,", ...",""))
        if(sum(colnames(D)%in%sID)<2) stop("Less than 2 samples from count table defined in sample meta")
    }
    D <- D[,sID,drop=F]
    return(D)
}
checkGeneName <- function(D,gID){
    if(sum(!gID%in%rownames(D))>0)
        stop(paste0("genes defined in count table (",length(gID),") are NOT available in effective length or TPM table (",nrow(D),")"))
    D <- D[gID,,drop=F]
    return(D)
}
getEffLength <- function(config,sID,gID,gInfo=NULL,gLength=NULL){
    D <- NULL
    if(!is.null(config$prj_effLength) && file.exists(config$prj_effLength)){
        message("reading effective length")
        #D <- read.table(config$prj_effLength,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
        D <- data.table::fread(config$prj_effLength,header=T,check.names=F)
        D <- data.frame(row.names=D[[1]],D[,-1],check.names=F)
        D <- checkSampleName(D,sID)
        D <- checkGeneName(D,gID)
        if(!is.null(gInfo) && !is.null(gLength)){
          selG <- rownames(D)[apply(D,1,median,na.rm=T)<gLength]
          if(length(selG)>0){
            message("\t\tThe nominal length of following ",length(selG)," genes will be used:")
            message("\t\t",paste(head(gInfo[selG,"Gene.Name"],5),collapse=", "),ifelse(length(selG)>5,", ...",""))
            D[selG,] <- matrix(gInfo[selG,"Length"],nrow=length(selG),ncol=ncol(D))
          }
        }
    }
    return(list(effLength=D))
}
getTPM <- function(config,sID,gID){
    D <- NULL
    if(!is.null(config$prj_TPM) && file.exists(config$prj_TPM)){
        message("reading sample TPM")
        #D <- read.table(config$prj_TPM,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
        D <- data.table::fread(config$prj_TPM,header=T,check.names=F)
        D <- data.frame(row.names=D[[1]],D[,-1],check.names=F)
        if(ncol(D)<2) stop("No sample detected!")
        D <- checkSampleName(D,sID)
        D <- checkGeneName(D,gID)
        D <- log2(config$count_prior+D)
    }
    return(list(logTPM=D))
}
getSeqQC <- function(config,sID){
    if(is.null(config$prj_seqQC)) return(NULL)
    message("reading sequence QC")
    #D <- read.table(config$prj_seqQC,sep="\t",header=T,as.is=T,check.names=F,row.names=1)
    D <- data.table::fread(config$prj_seqQC,header=T,check.names=F)
    D <- data.frame(row.names=D[[1]],D[,-1],check.names=F)
    if(ncol(D)<2) stop("No sample detected!")
    if(sum(!sID%in%rownames(D))>0)
        stop(paste0("samples (",paste(sID[!sID%in%rownames(D)],collapse=","),
                    ") defined in sample meta table are NOT available in sequence QC table"))
    D <- D[sID,,drop=F]
    return(list(seqQC=D))
}
filterGene <- function(config,D){
  if(length(config$rmGeneStartWith)>0 && !is.null(D$gInfo)){
    selRM <- rep(F,nrow(D$gInfo))
    for(one in config$rmGeneStartWith){
      selRM <- selRM | grepl(paste0("^",one),D$gInfo$Gene.Name)
    }
    if(sum(selRM)>0){
      message("The following genes will be filtered out:\n\t",
              paste(D$gInfo$Gene.Name[selRM],collapse=", "))
      selG <- D$gInfo$UniqueID[selRM]
      D$gInfo <- D$gInfo[!selRM,]
      D$counts <- D$counts[!rownames(D$counts)%in%selG,]
      if("logTPM"%in%names(D))
        D$logTPM <- D$logTPM[!rownames(D$logTPM)%in%selG,]
    }
  }
  return(D)
}
useAlias <- function(config,D){
    if(!is.null(config$sample_alias)){
        message("Applying alias")
        colnames(D$counts) <- rownames(D$meta) <- D$meta[,config$sample_alias]
        if(!is.null(D$effLength)) colnames(D$effLength) <- rownames(D$meta)
        if(!is.null(D$logTPM)) colnames(D$logTPM) <- rownames(D$meta)
        if(!is.null(D$seqQC)) rownames(D$seqQC) <- rownames(D$meta)
    }
    #print(names(D))
    return(D)
}
shortListMsg <- function(aVec,n=5){
    msg <- paste(head(aVec,n),collapse=", ")
    if(length(aVec)>n) msg <- paste(msg,", ... ")
    return(msg)
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
        ad_X <- rbind(ad_X,X[!ix,,drop=F])
        sizeF <- covariateRM_getSizeF(ad_X)
    }else if(grepl("combat",method,ignore.case=T)){
        ad_X <- covariateRM_ComBatRM(X[ix,],batchX[,1])
        ad_X <- rbind(ad_X,X[!ix,,drop=F])
        sizeF <- covariateRM_getSizeF(ad_X)
    }else{
        stop("unknown batch removal method!")
    }
    logTPM <- covariateRM_estTPM(ad_X,effeL[rownames(ad_X),],sizeF,prior)
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
setMetaFactor <- function(meta,strMetaFactor,addFactor=NULL,metaFactor=NULL){
    stopifnot(nrow(meta)>1)
    if(is.null(strMetaFactor))return(meta)
    ## retrieve the meta factor or initialize one ----
    meta <- metaFactor_checkNAempty(meta)
    if(is.null(metaFactor)){
        if(file.exists(strMetaFactor)){
            metaFactor <- yaml::read_yaml(strMetaFactor)
        }else{
            metaFactor <- metaFactor_addFactor(meta,strMetaFactor)
        }
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
            warning(paste0("The following information (",i,") defined in sample meta file are not defined in metafactor file:\n\t",
                           shortListMsg(oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]])))
            message("Please update metaFactor file to avoid this warning message")
            metaFactor[[i]] <- c(metaFactor[[i]],oneMetaUnique[!oneMetaUnique%in%metaFactor[[i]]])
        }
        if(sum(!metaFactor[[i]]%in%oneMetaUnique)>0){
            warning(paste0("The following information (",i,") defined in metafactor file are missing from sample file:\n\t",
                           shortListMsg(metaFactor[[i]][!metaFactor[[i]]%in%oneMetaUnique])))
            message("Please update metaFactor file to avoid this warning message")
            metaFactor[[i]] <- metaFactor[[i]][metaFactor[[i]]%in%oneMetaUnique]
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
    return(list(meta=meta,metaFactor=metaFactor))
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
## meta 'naive' ----
correctNaive <- function(config){
  for(one in c("sample_meta","sample_factor","comparison_file")){
    if(!file.exists(config[[one]])) next()
    a <- readLines(config[[one]])
    if(sum(grepl("ï",a,fixed=T))>0){
      message("'Naïve' detected in ",one,", replacing it")
      file.copy(config[[one]],paste0(config[[one]],".bak"))
      writeLines(gsub("ï","i",a,fixed=T),config[[one]])
    }
  }
}
## EAqc functions ------
source("PC_Covariates.R")
suppressPackageStartupMessages({
  require(ggplot2)
  require(reshape2)}) 
plotAlignQC <- function(estT,strPDF,estC=NULL,qc=NULL,prioQC=NULL,topN=c(1,10,30),gInfo=NULL,replot=F){#,50,100
    if(file.exists(strPDF) && !replot) return()
    message("plotting sequencing QC @",strPDF)
    if(!is.null(gInfo) & sum(rownames(gInfo)!=gInfo$Gene.Name,na.rm=T)>0){
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
plotMetaCor <- function(meta,strPDF){
    meta <- sapply(meta,function(x){
        if(length(unique(x))==1 || length(unique(x))==length(x)) return(NULL)
        return(x)
    })
    meta <- as.data.frame(meta[!sapply(meta,is.null)],stringsAsFactors = FALSE)
    pdf(strPDF)
    D <- plotEachCor(meta,meta)
    print(ggplot(D,aes(compGrp,CovGrp))+
              geom_tile(aes(fill=-log10(pvalue)))+
              xlab("Sample meta")+ylab("Sample meta")+ggtitle("log10 P-value")+
              scale_fill_gradient("-log10(pvalue) \n",low="#fee0d2",high="#99000d")+
              theme_minimal()+
              theme(axis.text.x=element_text(angle=90)))
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
    if(!is.null(selCov)) meta <- meta[,selCov,drop=F]
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
            logTPM <- suppressMessages(covariateRM(estC,effL,batchX=batchX,method=config$covariates_method,#'limma',
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
setSampleGroup <- function(config,meta){
    if(length(config$sample_group)>0)
        return(config)
    metaType <- sapply(meta,is.factor)
    for(one in names(metaType)){
        if(metaType[one] && nlevels(meta[[one]])>1 && nlevels(meta[[one]])<nrow(meta)){
            config$sample_group <- one
            return(config)
        }
    }
    config$sample_group <- names(metaType)[metaType][1]
    return(config)
}
finishQC <- function(strMsg){
    message("==========================================")
    system(paste0("cp ",strMsg$output,"/",strMsg$prj_name,"* ",strMsg$QuickOmics_test_folder))
    message(paste0("\n-----> Please visit: ",strMsg$QuickOmics_test_link,strMsg$prj_name))
    message("----->'EArun' can be used to obtain the QuickOmics object after necessary 'covariates_adjust' is set and comparison definition file is filled:")
    message("\t\t\t",strMsg$comparison_file)
    message("\t\tEArun ",strMsg$output,"/config.yml\n\n")
    message("-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.")
}

## EAsplit functions -------------
# EAdata: meta, counts, gInfo, effLength, logTPM, seqQC
splitUpdatePrj <- function(config,one,strOut){
    config$prj_name <- paste(config$prj_name,one,sep="_")
    config$prj_title <- paste0("(",one,") ",config$prj_title)
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
splitSaveFactor <- function(strSrc,strDest,strMeta,sep=","){
    if(file.exists(strSrc)){
        #meta <- read.table(strMeta,sep=sep,header=T,check.names=F,as.is=T)
        meta <- data.table::fread(strMeta,header=T,check.names=F)
        meta <- data.frame(row.names=meta[[1]],meta[,-1],check.names=F)
        metaF <- yaml::read_yaml(strSrc)
        conn <- file(strDest,"w")
        for(one in names(metaF)){
            res <- metaF[[one]][metaF[[one]]%in%meta[,one]]
            cat(one,": ['",paste(res,collapse="','"),"']\n",sep="",file=conn)
        }
        close(conn)
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
    config$prj_effLength <- splitSaveData(EAdata$effLength,
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
                                            config$sample_meta)
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
source("QuickOmics_DEG.R")

checkShinyTestSetting <- function(config){
  if(is.null(config$QuickOmics_test_folder) || 
     is.null(config$QuickOmics_test_link)){
    message("#####################")
    warning(paste("#################\nThe ShinyOne is not set up (No finishing web link)!\nPlease contact",
                  config$powerby,"\n#################\n"))
    return(F)
  }
  return(T)
}
saveCountsAlias <- function(config,estC){
    saveRDS(estC,file=paste0(config$output,"/",config$prj_name,"_estCount.rds"))
}
plotCovBio <- function(config,meta,comp_info){
  selCov <- unique(c(config$covariates_check,config$covariates_adjust))
  selCom <- unique(comp_info$Group_name)
  selCov <- selCov[!selCov%in%selCom]
  if(is.null(selCom) || is.null(selCov)) return()
  selAll <- c(selCov,selCom)
  meta <- meta[,selAll]

  ## change the Well_Row from charactor to numeric
  oneMeta <- "Well_Row"
  if(oneMeta %in% colnames(meta)) meta[,oneMeta] <- as.numeric(as.factor(meta[,oneMeta]))
  # remove meta with only one unique value
  meta <- meta[,apply(meta,2,function(x)return(length(unique(x))>1)),drop=F]
  selCov <- selCov[selCov%in%colnames(meta)]
  for(one in selCom) meta[,one] <- as.factor(meta[,one])
  
  message("Plotting correlation between covarites and comparison groups:\n\t",
          paste(selCov,collapse=", ")," v.s. ",paste(selCom,collapse=","))
  pdf(paste0(config$output,"/cov.vs.bio.pdf"))
  D <- plotEachCor(meta[,selCom,drop=F],
                   meta[,selCov,drop=F])
  print(ggplot(D,aes(compGrp,CovGrp))+
          geom_tile(aes(fill=-log10(pvalue)))+
          xlab("Comparison groups")+ylab("Covariate groups")+
          scale_fill_gradient("-log10(pvalue) \n",low="#fee0d2",high="#99000d")+
          theme_minimal()+
          theme(axis.text.x=element_text(angle=90)))
  a <- dev.off()
  data.table::fwrite(D,paste0(config$output,"/cov.vs.bio.csv"))
}
plotEachCor <- function(X,Y){
  corD <- NULL
  for(x in colnames(X)){
    for(y in colnames(Y)){
      message("\t",x," .vs. ",y)
      D <- cbind(X[,x,drop=F],Y[,y,drop=F]) %>% na.omit()
      oneCor <- data.frame(compGrp=x,CovGrp=y,method="",estimate=0,pvalue=1,stat=0,sampleN=nrow(D))
      if(x==y){
          corD <- rbind(corD,data.frame(compGrp=x,CovGrp=y,method="",estimate=1,pvalue=1e-10,stat=100,sampleN=nrow(D)))
          next
      }
      if(nrow(D)<3){
        corD <- rbind(corD,oneCor)
        next
      }else if(length(unique(D[,x]))<2 || length(unique(D[,y]))<2){
        message("\t\tSkip! only 1 unique value in ",x," or ",y," after NA removed.")
        corD <- rbind(corD,oneCor)
        next
      }else{
        oldNames <- colnames(D)
        colnames(D) <- names(oldNames) <- gsub("^[0-9]","x",gsub(" |[[:punct:]]","_",oldNames))
        x <- colnames(D)[1]
        y <- colnames(D)[2]
        if(is.numeric(D[,x]) && is.numeric(D[,y])){
          a <- cor.test(D[,x],D[,y],method="spearman")
          oneCor$method <- "spearman"
          oneCor$estimate <- as.vector(a$estimate)
          oneCor$pvalue <- as.vector(a$p.value)
          oneCor$stat <- as.vector(a$statistic)
          print(ggplot(D,aes_string(x,y))+
                  geom_point()+
                  geom_smooth()+
                  xlab(oldNames[x])+ylab(oldNames[y])+
                  theme_light())
        }else if(is.numeric(D[,x]) && (is.factor(D[,y]) || is.character(D[,y]))){
          a <- anova(lm(as.formula(paste(c(x,y),collapse="~")),data=D))
          oneCor$method <- "Anova"
          oneCor$pvalue <- as.vector(a$'Pr(>F)'[1])
          oneCor$stat <- as.vector(a$'F value'[1])
          suppressMessages(suppressWarnings({
          print(ggplot(D,aes_string(y,x))+
                    geom_violin(fill="yellow",color="gray",alpha=0.7) +
                    geom_jitter(width = 0.15,size=1,alpha=0.5,color="darkblue")+
                    xlab(oldNames[y])+ylab(oldNames[x])+
                    theme_light()+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
          }))
        }else if(is.numeric(D[,y]) && (is.factor(D[,x]) || is.character(D[,x]))){
          a <- anova(lm(as.formula(paste(c(y,x),collapse="~")),data=D))
          oneCor$method <- "Anova"
          oneCor$pvalue <- as.vector(a$'Pr(>F)'[1])
          oneCor$stat <- as.vector(a$'F value'[1])
          suppressMessages(suppressWarnings({
              print(ggplot(D,aes_string(x,y))+
                        geom_violin(fill="yellow",color="gray",alpha=0.7) +
                        geom_jitter(width = 0.15,size=1,alpha=0.5,color="darkblue")+
                        xlab(oldNames[x])+ylab(oldNames[y])+
                        theme_light()+
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
          }))
        }else if((is.factor(D[,x]) || is.character(D[,x])) && (is.factor(D[,y]) || is.character(D[,y]))){
          suppressMessages(suppressWarnings({
              require(gridExtra)
              require(grid)
          }))
          a <- suppressWarnings(chisq.test(D[,x],D[,y]))
          oneCor$method <- "Chisq"
          oneCor$pvalue <- as.vector(a$p.value)
          oneCor$stat <- as.vector(a$statistic)
          chiTable <- table(D)
          if(sum(nchar(rownames(chiTable)))<sum(nchar(colnames(chiTable)))) chiTable <- t(chiTable)
          labs <- names(dimnames(chiTable))
          
          legend_grob <- NULL
          rowCharN <- max(nchar(rownames(chiTable)))
          colCharN <- sum(nchar(colnames(chiTable)))
          charN_cutoff <- 45
          if((rowCharN+colCharN)>charN_cutoff){
              rLab <- setNames(paste0("R",1:nrow(chiTable)),rownames(chiTable))
              cLab <- setNames(paste0("C",1:ncol(chiTable)),colnames(chiTable))
              rowCharN1 <- max(nchar(rLab))
              colCharN1 <- sum(nchar(cLab))
              rLegend <- get_chiTable_legend(labs[1],rLab) 
              cLegend <- get_chiTable_legend(labs[2],cLab) 
              if((rowCharN1+colCharN)<charN_cutoff){
                  rownames(chiTable) <- rLab
                  legend_grob <- rLegend
              }else if((rowCharN+colCharN1)<charN_cutoff){
                  colnames(chiTable) <- cLab
                  legend_grob <- cLegend
              }else{
                  dimnames(chiTable) <- list(rLab,cLab)
                  legend_grob <- arrangeGrob(rLegend,cLegend,nrow=2)
              }
          }
          if(is.null(legend_grob)){
              grid.arrange(tableGrob(chiTable),
                           left=labs[1],
                           top=labs[2],
                           padding=unit(2, "line"))
          }else{
              grid.arrange(arrangeGrob(tableGrob(chiTable),
                                       left=labs[1],
                                       top=labs[2]),
                           legend_grob,
                           nrow=2,
                           heights=c(max(1,nrow(chiTable)*0.1),1),
                           padding=unit(2, "line"))
          }
        }
      }
      corD <- rbind(corD,oneCor)
    }
  }
  return(corD)
}
get_chiTable_legend <- function(title_text,labs,width=80){
    labs <- paste0(labs,": ",names(labs))
    nCol <- floor(width/max(nchar(labs)))
    if(nCol==1){
        labs <- list(labs)
    }else{
        labs <- split(labs,cut(seq_along(labs),
                               breaks = nCol,
                               labels = FALSE))
    }
    legend_labs <- lapply(1:length(labs),function(i){
        lab_title <- ""
        if(i==1) lab_title <- title_text 
        textGrob(paste(c(lab_title,labs[[i]]),
                       collapse="\n"),
                 x=0.05,y=1,just=c("left","top"),
                 gp=gpar(fontsize=10))
    })
    return(arrangeGrob(grobs=legend_labs,ncol=length(legend_labs)))
}
lowCountFiltering <- function(config,D,checkComp=T){
    if(is.null(config$minCounts)) return(D)
    message("Checking sample counts")
    sel <- apply(D$counts,2,sum) > config$minCounts
    if(sum(!sel)>0){
        message("The following samples are removed from DE analysis due to count minimal (",config$minCounts,") set in config (minCounts):")
        message(paste(c(head(colnames(D$counts)[!sel]),"..."),collapse=", "))
        message("===== Total of ",sum(sel)," samples after filtering")
        D$counts <- D$counts[,sel]
        D$meta <- D$meta[sel,]
        metaInfo <- setMetaFactor(D$meta,config$sample_factor,metaFactor=D$metaFactor)
        D$meta <- metaInfo$meta
        D$metaFactor <- metaInfo$metaFactor
        if(!is.null(D$effLength)) D$effLength <- D$effLength[,sel]
        if(!is.null(D$logTPM)) D$logTPM <- D$logTPM[,sel]
        if(!is.null(D$seqQC)) D$seqQC <- D$seqQC[sel,]
        # check the comparison informatio again after filtering
        if(checkComp)
            D$comp_info <- checkComparisonInfo(read_file(config$comparison_file,T),
                                               D$meta,config$comparison_file)
    }
    return(D)
}

comparisonAnalysis <- function(config,estC,meta,comp_info){
    message("====== Starting DE analysis ...")
    saveCountsAlias(config,estC)
    ## comparison -----------
    if((!is.null(config[['qsub']]) && config[['qsub']]) || (!is.null(config$parallel) && config$parallel=="sge")){
       message("sge DEG process ...")
        source(paste0(config$srcDir,"/qsubDEG.R"))
        return(qsubDEG(estC,meta,comp_info,config$output,config$srcDir,
                       core=config$core,qsubTime=config$qsubTime))
    }else if(!is.null(config$parallel) && config$parallel=="slurm"){
        message("slurm DEG process ...")
        source(paste0(config$srcDir,"/sbatchDEG.R"))
        return(sbatchDEG(estC,meta,comp_info,config$output,config$srcDir,
                     core=config$core))
    }else{
        message("serial DEG process ...")
        return(Batch_DEG(estC,meta,comp_info,core=config$core))
    }
}
plotDEG_MA <- function(deg,logTPM,meta,comp_info,config){
  suppressMessages(suppressWarnings(require(ggrastr)))
  pdf(paste0(config$output,"/deg.MA.pdf"))
  for(one in names(deg)){
    if(one %in% rownames(comp_info)){
      testID <- rownames(meta)[meta[,comp_info[one,"Group_name"]]==comp_info[one,"Group_test"]]
      ctrlID <- rownames(meta)[meta[,comp_info[one,"Group_name"]]==comp_info[one,"Group_ctrl"]]
      testExp <- apply(logTPM[,testID],1,mean)
      ctrlExp <- apply(logTPM[,ctrlID],1,mean)
      gID <- intersect(rownames(logTPM),rownames(deg[[one]]$DEG))
      X <- rbind(data.frame(gID=gID,
                            Exp=apply(cbind(testExp[gID],testExp[gID]),1,max),
                            log2FC=testExp[gID]-ctrlExp[gID],
                            method="from Exp"),
                 data.frame(gID=gID,
                            Exp=apply(cbind(testExp[gID],testExp[gID]),1,max),
                            log2FC=deg[[one]]$DEG[gID,grep("log2FoldChange",colnames(deg[[one]]$DEG))],
                            method="from DEG"))
      tryCatch(
        expr={
          p <- ggplot(X,aes(Exp,log2FC,color=method))+
            rasterise(geom_point(), dpi=100,dev="ragg")+
            scale_color_manual(values=c('from Exp'="#bdbdbd80",'from DEG'="#de2d26e0"))+
            ggtitle(one)+
            xlab("Expression (log2)")+
            ylab("log2FC")+
            theme_light()
          print(p)
        },
        error=function(e){
          p <- ggplot(X,aes(Exp,log2FC,color=method))+
            geom_point()+
            scale_color_manual(values=c('from Exp'="#bdbdbd80",'from DEG'="#de2d26e0"))+
            ggtitle(one)+
            xlab("Expression (log2)")+
            ylab("log2FC")+
            theme_light()
          print(p)
        }
      )
    }
  }
  a <- dev.off()
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
    if(is.null(X)){
        network <- data.frame(from=character(),to=character(),cor=numeric(),p=numeric(),direction=as.integer())
        save(network,file=paste0(config$output,"/",config$prj_name,"_network.RData"))
        return()
    }
    message("Obtaining networks ...")
	suppressMessages(suppressWarnings({
		require(Hmisc)
		require(tibble)
		}))
    cor_cutoff <- config$gene_network_cor_cutoff
    p_cutoff <- config$gene_network_p_cutoff
    variableN <- config$gene_network_high_variable_N
    edge_max <- as.numeric(config$gene_network_max_edge)
    edge_min <- as.numeric(config$gene_network_min_edge)

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
    if(length(config$sample_group)>0 && !"group" %in% colnames(meta)){
        meta <- cbind(meta,group=meta[,config$sample_group])
    }
    if(!"group" %in% colnames(meta)){
        comp <- read_file(config$comparison_file,T)
        meta <- setMetaFactor(meta,config$sample_factor,unique(comp$Group_name),metaFactor=D$metaFactor)$meta
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
        if(ncol(res)<5){
          res <- cbind(res,data.frame(row.names=rownames(res),
                                      rep(list(rep(NA,nrow(res))),5-ncol(res))))
        }
        colnames(res)[1:5] <- c("logFC","P.Value","Adj.P.Value","S.Value","Adj.S.Value")
        res <- res[,1:5]
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
  colnames(meta) <- gsub("^sampleid$","sample__id",
                         gsub("^Order$","Order__",
                              gsub("^ComparePairs$","Compare__Pairs",colnames(meta))))

  if(nrow(meta)>=length(comNames)){
    MetaData <- list(sampleid=rownames(meta),
                     group=meta$group,
                     Order=unique(meta$group),
                     ComparePairs=unlist(ifelse(is.null(comNames),"",list(comNames))))
  }else{
    MetaData <- list(sampleid=rownames(meta),
                     group=meta$group,
                     Order=unique(meta$group),
                     ComparePairs="")
  }
  MetaData <- as.data.frame(lapply(MetaData,'length<-',max(sapply(MetaData,length))),stringsAsFactors=F)
  MetaData$ComparePairs[is.na(MetaData$ComparePairs)] <- ""
  meta <- meta[,-which(colnames(meta)=="group"),drop=F]
  MetaData <- cbind(MetaData,meta)
  return(MetaData)
}
saveQuickOmics <- function(config,EAdata,DEGs=NULL){
    message("saving QuickOmics object ...")
    EAdata$meta <- updateMeta(config,EAdata$meta)
    comp_info <- read_file(config$comparison_file,T)
    
    message("\tFormating the expression data")
    data_wide <- EAdata$logTPM
    data_long <- melt(as.matrix(EAdata$logTPM))
    colnames(data_long) <- c("UniqueID","sampleid","expr")
    config$sample_group <- ifelse(length(config$sample_group)<1,"group",config$sample_group)
    data_long <- cbind(data_long,group=EAdata$meta[data_long$sampleid,config$sample_group])
    
    if(!is.null(DEGs)){
        message("\tFormating the DEG results")
        compRes <- formatQuickOmicsResult(DEGs,EAdata$logTPM,EAdata$meta[,"group"],EAdata$gInfo)
        data_results <- compRes$Dw
        results_long <- compRes$Dl
    }else{
        data_results <- NULL
        comp_info <- NULL
        results_long <- data.frame(UniqueID=character(),test=factor(),logFC=numeric(),P.Value=numeric(),
                                   Adj.P.Value=numeric(),S.Value=numeric(),Adj.S.Value=numeric())#column type matters in Quickomics
    }
    
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
                         Path=config$output,
                         CovariantMethod=ifelse(length(config$covariates_adjust)>0,
                                                config$covariates_method,"")),
              file=paste0(config$output,"/",config$prj_name,".csv"),
              row.names=F)
}
finishRun <- function(strMsg){
    message("=================================================\nResults are saved in ",strMsg$output)
    if(!checkShinyTestSetting(sysConfig)) return()
    system(paste0("cp ",strMsg$output,"/",strMsg$prj_name,"* ",strMsg$QuickOmics_test_folder))
    message(paste0("\n-----> Please visit: ",strMsg$QuickOmics_test_link,strMsg$prj_name))
    message("Please carefully review the results before publishing:")
    message("----->'EApub' can be used to publish the project into ShinyOne project manager: ",
            sapply(strsplit(strMsg$shinyApp,"\\/"),
                   function(x)return(paste(x[1:grep("shiny",x)],collapse="/"))))
    
}

## EApub functions --------
checkShinySetting <- function(config){
  if(is.null(config$QuickOmics_test_folder) || 
     is.null(config$QuickOmics_test_link) || 
     is.null(config$QuickOmics_publish_folder) || 
     is.null(config$QuiclOmics_publish_link) || 
     is.null(config$shinyApp) || 
     is.null(config$shinyAppKey))
    stop(paste("The ShinyOne is not set up!\nPlease contact",config$powerby))
}
checkPubFiles <- function(config){
    if(!file.exists(paste0(config$output,"/",config$prj_name,".RData")))
        stop("Missing project Quickomics object")
    if(!file.exists(paste0(config$output,"/",config$prj_name,"_network.RData")))
        stop("Missing project Quickomics network object")
    if(!file.exists(paste0(config$output,"/",config$prj_name,".csv")))
        stop("Missing project discription csv file")
}
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
    shinyOneData[["Data_Cleaning"]] <- paste("Data processing folder:",config$output)
    return(shinyOneData)
}
pubShinyOne <- function(config){
    checkShinySetting(config)
    checkPubFiles(config)
    strF <- paste0(config$QuickOmics_publish_folder,config$prj_name,".RData")
    if(file.exists(strF)){
        message("project files exists in ",config$QuickOmics_publish_folder)
        stop("The project already exists in ShinyOne!\nPlease remove the record and associated files or change prj_name and rerun EApub!")
    }
    message("preparing information for ShinyOne")
    shinyOneData <- getShinyOneInfo(config)
    shinyOneCMD <- paste0("curl -s -k -X POST -d 'data={",
                          paste(paste0('"',names(shinyOneData),'": "',shinyOneData,'"'),collapse = ", "),
                          "}' '",config$shinyApp,"api_add_project.php?api_key=",config$shinyAppKey,"'")

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
                    "\nPlease contact ",config$powerby))
    }
    system(paste0("cp ",config$output,"/",config$prj_name,"* ",config$QuickOmics_publish_folder))
    return(shinyMsg$ID)
}
finishShinyOne <- function(shinyMsg){
    message("=================================================\nShinyOne access: ",
            shinyMsg$shinyApp,"app_project_review.php?ID=",shinyMsg$ID)
    
}

## CellMap functions -----
cellmap_avaid_profiles <- function(sysConfig){
    cm_profiles <- NULL
    if(is.null(sysConfig)) return(cm_profiles)
    if(!is.null(sysConfig$CellMapExe) && file.exists(sysConfig$CellMapExe)){
        res <- system(sysConfig$CellMapExe,intern=T,ignore.stderr=T)
        iStart <- grep("-p PROFILE",res)+2
        iEnd <- seq_along(res)[nchar(res)==0]
        iEnd <- iEnd[iEnd>iStart][1]-1
        if(iEnd>iStart){
            message("*** CellMap is avaialble with the following profiles: ***")
            for(i in iStart:iEnd){
                pf <- sapply(strsplit(res[i],":"),head,1)
                if(nchar(pf)<15){
                    message(res[i])
                    cm_profiles <- c(cm_profiles,trimws(pf))
                }
            }
        }
    }
    return(cm_profiles)
}
#CellMap_profile
cellmap_run <- function(logTPM,config,sysConfig){
    if(!is.null(config$CellMap_profile) && !is.null(sysConfig$CellMapExe) && file.exists(sysConfig$CellMapExe)){
        message("*** CellMap is running on the cluster, please check the CellMap sub-folder in the output folder ***")
        TPM <- 2^logTPM-config$count_prior
        strPrefix <- file.path(config$output,"cellmap",config$prj_name)
        dir.create(dirname(strPrefix),showWarnings=F)
        strTPM <- paste0(strPrefix,"_TPM.tsv")
        write.table(TPM,file=strTPM,sep="\t")
        strCMD <- paste(sysConfig$CellMapExe,"-b",strTPM,"-p",config$CellMap_profile,"-o",strPrefix)
        if(!is.null(config$CellMap_rm_ct) && length(config$CellMap_rm_ct)>0)
            strCMD <- paste0(strCMD," -d ",paste(config$CellMap_rm_ct,collapse=","))
        if(!is.null(config$CellMap_species) && length(config$CellMap_species)>0)
            strCMD <- paste(strCMD,"--species",config$CellMap_species)
        system(strCMD)
    }
}
## others ----------
initialMsg(dirname(getwd()))
