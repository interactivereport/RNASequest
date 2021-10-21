#!/usr/bin/env Rscript
rm(list=ls())
.libPaths(grep("home",.libPaths(),invert=T,value=T))
args = commandArgs(trailingOnly=T)
## please note: the command section started with EApub is a indicator this section of code will be published by EApub
## if/else right after indicating which section of code will be included,
## thus such section will only contain one pair of if/else
## init -----
if(length(args)<2){
    message("An example config can be found: /camhpc/ngs/projects/TST11589/dnanexus/20210426220540_Zhengyu.Ouyang/config.yml")
    message("'EAinit' can be used to create a config file for a RNAseq project")
    stop("config yaml file is required!")
}
message("Loading resources ...")
config <- sapply(yaml::read_yaml(args[2]),unlist)
sys_config <- yaml::read_yaml(paste0(args[1],"sys.yml"))
source(paste0(args[1],"covariateRM.R"))
source(paste0(args[1],"QuickOmics_DEG.R"))
source(paste0(args[1],"formatQuickOmicsResult.R"))
source(paste0(args[1],"getNetwork.R"))
source(paste0(args[1],"Hmisc.rcorr.R"))
source(paste0(args[1],"readData.R"))
source(paste0(args[1],"infoCheck.R"))
source(paste0(args[1],"qsubDEG.R"))
source(paste0(args[1],"metaFactor.R"))

checkConfig(config)

## EApub: read and check the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,check.names=F,as.is=T)
checkConsistConfigMeta(config,meta)
rownames(meta) <- meta[,config$sample_name]

## EApub: read and check the comparison file ------
message("====== reading comparison information ...")
comp_info <- checkComparisonInfo(read_file(config$comparison_file,T),
                                 meta,config$comparison_file)
## EApub: set meta factors ------
meta <- metaFactor(meta,config$sample_factor,unique(comp_info$Group_name))

## EApub: use first group name in comparison file for group information if group is not defined in meta -----
if(!'group'%in%colnames(meta)){
    if(is.null(config$sample_group) || length(config$sample_group)==0){
        config$sample_group <- comp_info[1,"Group_name"]
    }
    meta <- cbind(group=apply(meta[,config$sample_group,drop=F],1,function(x)return(paste(x,sep="."))),meta)
}

## EApub if: gene definition file ---------
message("====== reading gene annotation ...")
if(!is.null(config$gene_annotation)){
    gInfo <- read.csv(config$gene_annotation,row.names=1,as.is=T)
}else{
    gInfo <- data.frame(id=rownames(estCount),
                        UniqueID=rownames(estCount),
                        Gene.Name=rownames(estCount),
                        Biotype='unknown')
}

## EApub if: read the gene quantification input ----
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
## EApub: process gene quantification -----
checkSampleName(rownames(meta),colnames(estCount))
if(!is.null(estCount)) estCount <- estCount[,rownames(meta)]
if(!is.null(effeL)) effeL <- effeL[,rownames(meta)]
estCount <- estCount[apply(estCount,1,function(x)return(sum(x>=config$min_count)))>=config$min_sample,]

## EApub: covariates removal ----
message("====== adjusting covariates for visualization ...")
if(!is.null(estCount) && !is.null(effeL)){
    if(is.null(config$covariates_adjust)){
        message("'covariates_adjust' is not set in the config file, no covariate adjust")
        batchX <- NULL
        yaxisLab <- paste0("log2(TPM+",config$count_prior,")")
    }
    else{
        batchX <- meta[,config$covariates_adjust,drop=F]
        yaxisLab <- paste0("log2(estTPM+",config$count_prior,")")
    }
    logTPM <- covariateRM(estCount,effeL,batchX=batchX,method='limma',
                          prior=config$count_prior)
}
if(is.null(logTPM)){
    stop("Gene quantification is missing, please provide either raw counts with effective length or TPM")
}
logTPM <- logTPM[rownames(estCount),rownames(meta)]
## EApub: sample alias ------
if(!is.null(config$sample_alias)){
    rownames(meta) <- colnames(logTPM) <- colnames(estCount) <- meta[,config$sample_alias]
}
saveRDS(estCount,file=paste0(config$output,"/",config$prj_name,"_estCount.rds"))
## EApub else: comparison -----------
message("====== Starting DEG analyses ...")
if(!is.null(config$qsub) && config$qsub){
    DEGs <- qsubDEG(estCount,meta,comp_info,config$output,args[1],core=config$core)
}else{
    DEGs <- Batch_DEG(estCount, meta, comp_info,core=config$core)
}
message("Formating the DEG results")
compRes <- formatQuickOmicsResult(DEGs,logTPM,meta[,"group"],gInfo)
data_results <- compRes$Dw
results_long <- compRes$Dl
## EApub: produce the network for quickOmics ----------
message("====== gene network generation ...")
network <- getNetwork(logTPM,
                      cor_cutoff=config$gene_network_cor_cutoff,
                      p_cutoff=config$gene_network_p_cutoff,
                      variableN=config$gene_network_high_variable_N,
                      edge_max=as.numeric(config$gene_network_max_edge),
                      edge_min=as.numeric(config$gene_network_min_edge),
                      core=config$core)
save(network,file=paste0(config$output,"/",config$prj_name,"_network.RData"))
## EApub: save the R object for quickOmics--------
message("====== saving QuickOmics object ...")
data_wide <- logTPM
data_long <- melt(as.matrix(logTPM))
colnames(data_long) <- c("UniqueID","sampleid","expr")
data_long <- cbind(data_long,group=meta[data_long$sampleid,config$sample_group])


MetaData <- formatQuickOmicsMeta(meta,names(DEGs))
ProteinGeneName <- gInfo[rownames(logTPM),]
save(data_results,results_long,
     data_wide,data_long,
     MetaData,ProteinGeneName,
     comp_info,
     yaxisLab,
     file=paste0(config$output,"/",config$prj_name,".RData"))
## save the project csv file -------
write.csv(data.frame(Name=ifelse(is.null(config$prj_title),
                                 config$prj_name,
                                 paste(config$prj_name,config$prj_title,sep=": ")),
                     ShortName=config$prj_name,
                     ProjectID=config$prj_name,
                     Species=config$species, 
                     ExpressionUnit=yaxisLab,
                     Path=config$output),
          file=paste0(config$output,"/",config$prj_name,".csv"),
          row.names=F)

## finishing -----
#system(paste0("cp ",config$output,"/",config$prj_name,"* "))
message("=================================================\nResults are saved in ",config$output)
system(paste0("cp ",config$output,"/",config$prj_name,"* ",sys_config$QuickOmics_test_folder))
#message("-----> to check the results on QuickOmics, execute the following command on ngs server")
#message("\t\tcp ",config$output,"/",config$prj_name,"* /srv/shiny-server/Quickomics/unlisted/.")
#message("\t\t\t And then Please visit: http://ngs.biogen.com:3838/Quickomics/?unlisted=",config$prj_name)
message(paste0("\n-----> Please visit: ",sys_config$QuickOmics_test_link,config$prj_name))

message("\nPowered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")

sink(paste0(config$output,"/session.EArun"))
sessionInfo()
sink()
