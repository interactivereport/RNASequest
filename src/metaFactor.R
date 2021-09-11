
metaFactor <- function(meta,strMetaFactor,addFactor=NULL){
    if(is.null(strMetaFactor))return(meta)
    if(file.exists(strMetaFactor)){
        metaFactor <- yaml::read_yaml(strMetaFactor)
    }else{
        metaFactor <- list()
        for(i in colnames(meta)){
            if(grepl("URL",i)) next
            if(is.character(meta[,i]) || is.factor(meta))
                metaFactor[[i]] <- unique(meta[,i])
        }
        saveYaml(metaFactor,strMetaFactor)
    }
    if(sum(!names(metaFactor)%in%colnames(meta))>0){
        stop(paste0("unknown meta (",
                    paste(names(metaFactor)[!names(metaFactor)%in%colnames(meta)],collapse=","),
                    ") defined in meta factor file(",strMetaFactor,")"))
    }
    for(i in names(metaFactor)){
        if(sum(!unique(meta[,i])%in%metaFactor[[i]])>0){
            stop(paste0(paste(unique(meta[,i])[!unique(meta[,i])%in%metaFactor[[i]]],collapse=","),
                        " from ",i," are not defined in meta factor file"))
        }
        meta[,i] <- factor(meta[,i],levels = metaFactor[[i]])
    }
    if(!is.null(addFactor)){
        for(i in addFactor){
            if(!is.factor(meta[,i])){
                meta[,i] <- factor(meta[,i],levels=unique(meta[,i]))
            }
        }
    }
    return(meta)
}

saveYaml <- function(ymlist,strF){
    cat(paste(paste0(names(ymlist),": ['",
                     sapply(ymlist,paste,collapse="','"),"']"),
              collapse="\n"),
        "\n",sep="",file=strF)
}