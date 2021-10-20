
metaFactor <- function(meta,strMetaFactor,addFactor=NULL){
    if(is.null(strMetaFactor))return(meta)
    ## retrieve the meta factor or initialize one ----
    meta <- checkNAempty(meta)
    if(file.exists(strMetaFactor)){
        metaFactor <- yaml::read_yaml(strMetaFactor)
    }else{
        metaFactor <- addFactor(meta,strMetaFactor)
    }
    ## if any annotation removed from meta but existing in the meta factor ----
    if(sum(!names(metaFactor)%in%colnames(meta))>0){
        warning("The following removed from sampleMeta table, updating metaFactor file:")
        message(paste(names(metaFactor)[!names(metaFactor)%in%colnames(meta)],collapse=","))
        metaFactor <- metaFactor[names(metaFactor)%in%colnames(meta)]
        saveYaml(metaFactor,strMetaFactor)
    }
    ## add additional (new) meta annotation into the meta factor -------
    # since there are numerical annotation, the following mostly always will be executed
    if(sum(!colnames(meta)%in%names(metaFactor))>0)
        metaFactor <- addFactor(meta,strMetaFactor,metaFactor)

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

checkNAempty <- function(meta){
    for(i in colnames(meta)){
        if(grepl("URL",i)) next
        if(is.character(meta[,i]) || is.factor(meta)){
            meta[,i] <- as.character(meta[,i])
            meta[is.na(meta[,i])|nchar(meta[,i])==0,i] <- "NA"
        }
    }
    return(meta)
}

addFactor <- function(meta,strMetaFactor,metaFactor=list()){
    for(i in colnames(meta)){
        if(grepl("URL",i)) next
        if(is.character(meta[,i]) || is.factor(meta))
            metaFactor[[i]] <- unique(meta[,i])
    }
    saveYaml(metaFactor,strMetaFactor)
    return(metaFactor)
}

saveYaml <- function(ymlist,strF){
    cat(paste(paste0(names(ymlist),": ['",
                     sapply(ymlist,paste,collapse="','"),"']"),
              collapse="\n"),
        "\n",sep="",file=strF)
}