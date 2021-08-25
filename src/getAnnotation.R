getAnnotation <- function(strF,strPath=NULL){
    if(grepl("json$",strF)){
        gConfig <- rjson::fromJSON(file=strF)
        #gtfPath <- paste0(strPath,"/rnaseq/",
        #                  gConfig$global_params$reference$species,
        #                  "/",gConfig$global_params$reference$version,
        #                  "/",gConfig$global_params$reference$version,
        #                  ".transcript.gtf.gz")
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
    return(gConfig$global_params$internal_project_id)
}