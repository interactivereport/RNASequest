## check config file -----
checkConfig <- function(config){
    if(is.null(config$sample_name)) stop("sample_name needs to be defined in config.yml. Default is Sample_Name.")
    if(!file.exists(config$sample_meta)) stop(paste("The sample meta file, ",config$sample_meta,", does NOT exist!"))
    if(!file.exists(config$gene_annotation)) stop(paste("The gene annotation file, ",config$gene_annotation,", does NOT exist!"))
    if(!file.exists(config$comparison_file)) stop(paste("The comparison definition file, ",config$comparison_file,", does NOT exist!"))
}
## check consistency between config file and meta file -----
checkConsistConfigMeta <- function(config,meta){
    if(sum(config$sample_name==colnames(meta))!=1)
        stop(paste0("sample_name (",config$sample_name,") in config is NOT defined in the sample meta file."))
    if(!is.null(config$sample_alias)){
        if(sum(config$sample_alias%in%colnames(meta))==0)
            stop(paste("alias column",config$sample_alias,"is NOT defined in the sample meta file"))
        if(sum(duplicated(meta[,config$sample_alias]))>0)
            stop(paste("alias column",config$sample_alias,"contains duplicates in the sample meta file"))
    }
    if(!is.null(config$split_meta) && sum(config$split_meta==colnames(meta))!=1)
        stop(paste0("split_meta (",config$split_meta,") in config is NOT defined in the sample meta file."))
    
    if(!is.null(config$covariates_check) && sum(!config$covariates_check%in%colnames(meta))>0)
        stop(paste0("covariates_check variables defined in config (",
                    paste(config$covariates_check[!config$covariates_check%in%colnames(meta)],collapse=", "),
                   "), is NOT defined in the sample meta file"))
    if(!is.null(config$covariates_adjust) && sum(!config$covariates_adjust%in%colnames(meta))>0)
        stop(paste0("covariates_adjust variables defined in config (",
                    config$covariates_adjust[!config$covariates_adjust%in%colnames(meta)],
                    "), is NOT defined in the sample meta file"))

}
## check the same name consistency between meta and quantify ------
checkSampleName <- function(mID,qID){
    if(sum(!mID%in%qID)>0){
        stop(paste("sample names in the meta table do not match in the quantification table:\n",paste(mID[!mID%in%qID],collapse="\t")))
    }
    
}