#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=T)
## init -----
if(length(args)<2){
    message("An example config can be found: /camhpc/ngs/projects/TST11589/dnanexus/20210426220540_Zhengyu.Ouyang/config.yml")
    message("'EAinit' can be used to create a config file for a RNAseq project")
    stop("config yaml file is required!")
}
message("Loading resources ...")
config <- yaml::read_yaml(args[2])
#sys_config <- yaml::read_yaml(paste0(args[1],"sys.yml"))
source(paste0(args[1],"PC_Covariates.R"))
source(paste0(args[1],"covariateRM.R"))

system(paste0("rm -f ",config$output,"/covariatePCanalysis_*"))
## read the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,row.names=1,check.names=F,as.is=T)
if(!is.null(config$sample_name))
    rownames(meta) <- meta[,config$sample_name]
meta <- meta[,unique(c(config$covariates_check,config$covariates_adjust))]
## change the Well_Row from charactor to numeric
oneMeta <- "Well_Row"
if(oneMeta %in% colnames(meta)) meta[,oneMeta] <- as.numeric(as.factor(meta[,oneMeta]))
## read the gene quantification input ----
message("====== reading gene quantification ...")
estCount <- effeL <- logTPM <- yaxisLab <- NULL
if(!is.null(config$prj_path)){
    estCount <- read.table(paste0(config$prj_path,"/combine_rsem_outputs/genes.estcount_table.txt"),
                           header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    estCount <- estCount[!grepl("^ERCC",rownames(estCount)),]
    effeL <- read.table(paste0(config$prj_path,"/combine_rsem_outputs/genes.effective_length.txt"),
                        header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    effeL <- effeL[!grepl("^ERCC",rownames(effeL)),]
    colnames(estCount) <- sapply(strsplit(sapply(strsplit(colnames(estCount),"\\|"),
                                                 function(x)return(paste(head(x,-1),
                                                                         collapse="|"))),
                                          "_"),function(x)return(paste(x[-1],
                                                                       collapse="_")))
    colnames(effeL) <- sapply(strsplit(sapply(strsplit(colnames(effeL),"\\|"),
                                              function(x)return(paste(head(x,-1),
                                                                      collapse="|"))),
                                       "_"),function(x)return(paste(x[-1],
                                                                    collapse="_")))
}else{
    if(!is.null(config$exp_counts))
        estCount <- read.table(config$exp_counts,
                               header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    if(!is.null(config$exp_effective_length))
        effeL <- read.table(config$exp_effective_length,
                            header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    if(!is.null(config$exp_tpm)){
        logTPM <- log2(config$count_prior+
                           read.table(config$exp_tpm,header=T,row.names=1,
                                      sep="\t",check.names=F,as.is=T))
        yaxisLab <- paste0("log2(TPM+",config$count_prior,")")
    }
}
if(sum(!rownames(meta)%in%colnames(estCount))>0){
    message("The following samles from meta information is not in count matrix:")
    message(paste(rownames(meta)[!rownames(meta)%in%colnames(estCount)],collapse="\n"))
    stop("sample names in the meta table do not match sample name in the count table!")
}
if(!is.null(estCount)) estCount <- estCount[,rownames(meta)]
if(!is.null(effeL)) effeL <- effeL[,rownames(meta)]
estCount <- estCount[apply(estCount,1,function(x)return(sum(x>=config$min_count)))>=config$min_sample,]

## PCA QC analysis before covariates removal ---------
message("====== TPM estimation ...")
logTPM <- covariateRM(estCount,effeL,method=NULL)
res <- suppressMessages(suppressWarnings(
    Covariate_PC_Analysis(logTPM,meta,
                          out_prefix=paste0(config$output,"/covariatePCanalysis_noAdjust"),
                          PC_cutoff=config$covariates_check_PCcutoff,
                          FDR_cutoff=config$covariates_check_FDRcutoff,
                          N_col=config$covariates_check_plotNcol)))
message("============================================================")
message("-----> PC analysis is done with no covariate adjusted:\n\t",config$output,"/covariatePCanalysis_noAdjust...")
## covariates removal ----
if(is.null(config$covariates_adjust) || length(config$covariates_adjust)==0){
    warning("< covariates_adjust is NOT set in the config file, no covariate was adjusted! >")
}else{
    message("====== removing covariates for visualization ...")
    if(!is.null(estCount) && !is.null(effeL)){
        batchX <- meta[,config$covariates_adjust,drop=F]
        logTPM <- suppressMessages(covariateRM(estCount,effeL,batchX=batchX,method='limma',
                              prior=config$count_prio))
        #save(logTPM,meta,file="PCanalysis.rdata")
        res <- suppressMessages(suppressWarnings(
            Covariate_PC_Analysis(logTPM,meta,
                                  out_prefix=paste0(config$output,"covariatePCanalysis_Adjusted"),
                                  PC_cutoff=config$covariates_check_PCcutoff,
                                  FDR_cutoff=config$covariates_check_FDRcutoff,
                                  N_col=config$covariates_check_plotNcol)))
        message("-----> PC analysis is done after covariate removal:\n\t",config$output,"/covariatePCanalysis_Adjusted...")
        
    }else{
        warning("effective length file is NOT available, no covariate removing!")
    }
}

## finishing --------
message("==========================================")
message("----->'EArun' can be used to obtain the QuickOmics object after necessary 'covariates_adjust' is set and comparison definition file is filled:")
message("\t\t\t",config$comparison_file)
message("\t\tEArun ",config$output,"/config.yml")
message("Powered by the Computational Biology Group [zhengyu.ouyang@biogen.com]")

sink(paste0(config$output,"/session.EAqc"))
sessionInfo()
sink()




