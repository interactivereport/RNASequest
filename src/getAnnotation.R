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
                              "gtf.gz$",full.names=T)
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
    return(gInfo)
}

getSpecies <- function(strF){
    gConfig <- rjson::fromJSON(file=strF)
    return(gConfig$global_params$reference$species)
}
getProjectName <- function(strF){
    gConfig <- rjson::fromJSON(file=strF)
    return(gConfig$global_params$internal_project_id)
}