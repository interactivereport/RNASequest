#!/usr/bin/env Rscript
rm(list=ls())
.libPaths(grep("home",.libPaths(),invert=T,value=T))
args = commandArgs(trailingOnly=T)
## init -----
if(length(args)<2){
    message("An example config can be found: /camhpc/ngs/projects/TST11589/dnanexus/20210426220540_Zhengyu.Ouyang/config.yml")
    message("'EAinit' can be used to create a config file for a RNAseq project")
    args <- c("/home/zouyang/projects/quickOmics/src/","/camhpc/ngs/projects/TST11781/dnanexus/20210607022720_Maria.Zavodszky/EA20210902_0/Sc/config.yml")
    #stop("config yaml file is required!")
}
message("Loading resources ...")
config <- sapply(yaml::read_yaml(args[2]),unlist)
sys_config <- yaml::read_yaml(paste0(args[1],"sys.yml"))
source(paste0(args[1],"PC_Covariates.R"))
source(paste0(args[1],"covariateRM.R"))
source(paste0(args[1],"alignQC.R"))
source(paste0(args[1],"readData.R"))
source(paste0(args[1],"infoCheck.R"))
source(paste0(args[1],"metaFactor.R"))

checkConfig(config)
system(paste0("rm -f ",config$output,"/covariatePCanalysis_*"))
## read the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,check.names=F,as.is=T)#,row.names=1
checkConsistConfigMeta(config,meta)
rownames(meta) <- meta[,config$sample_name]
## set meta factors ------
meta <- metaFactor(meta,config$sample_factor)
## read the gene annotation -----
gInfo <- read.csv(config$gene_annotation,row.names=1,as.is=T)
## read the gene quantification input ----
message("====== reading gene quantification ...")
estCount <- effeL <- logTPM <- yaxisLab <- NULL
if(!is.null(config$prj_path)){
    estCount <- readData(paste0(config$prj_path,"/combine_rsem_outputs/genes.estcount_table.txt"),
                         rownames(meta),rownames(gInfo))
    effeL <- readData(paste0(config$prj_path,"/combine_rsem_outputs/genes.effective_length.txt"),
                      rownames(meta),rownames(gInfo))
    if(config$min_median_effective_length>0){
        message("\tfiltering by effective gene length>=",config$min_median_effective_length)
        effeL <- effeL[apply(effeL,1,median)>=config$min_median_effective_length,]
        estCount <- estCount[rownames(effeL),]
        message("\t\t",nrow(estCount)," genes")
    }
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
estCount <- estCount[apply(estCount,1,function(x)return(sum(x>=config$min_count)))>=config$min_sample,]

## read alignment QC and plot alignment QC if alias exists ----
qc <- readQC(paste0(config$prj_path,"/combine_rnaseqc/combined.metrics.tsv"),rownames(meta))
if(!is.null(config$sample_alias)){
    message("====== Plot alignment QC for alias ...")
    logTPM <- covariateRM(estCount,effeL,method=NULL,prior=config$count_prior)
    estT <- 2^logTPM-config$count_prior
    rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
    estC <- estCount
    rownames(estC) <- paste(rownames(estC),gInfo[rownames(estC),"Gene.Name"],sep="|")
    alignQC(estT,
            qc,
            paste0(config$output,"/alignQC.alias.pdf"),
            prioQC=sys_config$qc2meta,
            sIDalias=setNames(meta[,config$sample_alias],rownames(meta)),
            estC=estC)
}
## select covariates for analysis-------
selCov <- unique(c(config$covariates_check,config$covariates_adjust))
meta <- meta[,selCov]

## change the Well_Row from charactor to numeric
oneMeta <- "Well_Row"
if(oneMeta %in% colnames(meta)) meta[,oneMeta] <- as.numeric(as.factor(meta[,oneMeta]))

## PCA QC analysis before covariates removal ---------
message("====== TPM estimation ...")
logTPM <- covariateRM(estCount,effeL,method=NULL,prior=config$count_prior)
res <- suppressMessages(suppressWarnings(
    Covariate_PC_Analysis(logTPM,meta,
                          out_prefix=paste0(config$output,"/notAdjusted"),
                          PC_cutoff=config$covariates_check_PCcutoff,
                          FDR_cutoff=config$covariates_check_FDRcutoff,
                          N_col=config$covariates_check_plotNcol)))
message("============================================================")
message("-----> PC analysis is done with no covariate adjusted:\n\t",config$output,"/covariatePCanalysis_noAdjust...")
## PCA QC analysis after covariates removal ----
if(is.null(config$covariates_adjust) || length(config$covariates_adjust)==0){
    warning("< covariates_adjust is NOT set in the config file, no covariate was adjusted! >")
}else{
    message("====== removing covariates for visualization ...")
    if(!is.null(estCount) && !is.null(effeL)){
        batchX <- meta[,config$covariates_adjust,drop=F]
        logTPM <- suppressMessages(covariateRM(estCount,effeL,batchX=batchX,method='limma',
                              prior=config$count_prior))
        estT <- 2^logTPM-config$count_prior
        rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
        estC <- estCount
        rownames(estC) <- paste(rownames(estC),gInfo[rownames(estC),"Gene.Name"],sep="|")
        if(!is.null(config$sample_alias))
            alignQC(estT,
                    qc,
                    paste0(config$output,"/Adjusted.alignQC.pdf"),
                    prioQC=sys_config$qc2meta,
                    sIDalias=setNames(meta[,config$sample_alias],rownames(meta)),
                    estC=estC)
        else
            alignQC(estT,
                    qc,
                    paste0(config$output,"/Adjusted.alignQC.pdf"),
                    prioQC=sys_config$qc2meta,
                    estC=estC)
        res <- suppressMessages(suppressWarnings(
            Covariate_PC_Analysis(logTPM,meta,
                                  out_prefix=paste0(config$output,"/Adjusted"),
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
message("\t\tEArun ",config$output,"/config.yml\n\n")
message("-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.\n")

message("Powered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")

sink(paste0(config$output,"/session.EAqc"))
sessionInfo()
sink()
