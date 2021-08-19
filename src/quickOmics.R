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
config <- sapply(yaml::read_yaml(args[2]),unlist)
sys_config <- yaml::read_yaml(paste0(args[1],"sys.yml"))
source(paste0(args[1],"covariateRM.R"))
source(paste0(args[1],"QuickOmics_DEG.R"))
source(paste0(args[1],"formatQuickOmicsResult.R"))
source(paste0(args[1],"getNetwork.R"))
source(paste0(args[1],"Hmisc.rcorr.R"))
source(paste0(args[1],"readData.R"))
source(paste0(args[1],"infoCheck.R"))

checkConfig(config)
## check the comparison file, stop if empty ----------
comp_info <- read_file(config$comparison_file,T)
if(nrow(comp_info)==0){
    stop(paste0("Empty comparison definition file (",config$comparison_file,") is NOT allowed!"))
}
if(is.null(config$sample_group) || length(config$sample_group)==0){
    config$sample_group <- comp_info[1,"Group_name"]
}

# set the default value if those are empty
setDefault <- F
for(i in rownames(comp_info)){
    if(is.null(comp_info[i,"Group_name"]) || nchar(comp_info[i,"Group_name"])==0){
        stop(paste0("'Group_name' cannot be empty for comparison",i))
    }
    if(is.null(comp_info[i,"Model"]) || nchar(comp_info[i,"Model"])==0){
        comp_info[i,"Model"] <- comp_info[i,"Group_name"]
        setDefault <- T
    }
    if(!grepl(comp_info[i,"Group_name"],comp_info[i,"Model"])){
        stop(paste("The DEG 'Model' didn't include 'Group_name' in the comparison file for",i))
    }
    if(is.null(comp_info[i,"Shrink_logFC"]) || nchar(comp_info[i,"Shrink_logFC"])==0){
        comp_info[i,"Shrink_logFC"] <- "Yes"
        setDefault <- T
    }
    if(is.null(comp_info[i,"LFC_cutoff"]) || nchar(comp_info[i,"LFC_cutoff"])==0){
        comp_info[i,"LFC_cutoff"] <- 0
        setDefault <- T
    }
    if(!comp_info[i,"Analysis_method"]%in%c("DESeq2","limma")){
        stop(paste(comp_info[i,"Analysis_method"],"for comparison",i,
                   "is NOT a valide 'Analysis_method' (DESeq2 or limma) in comparison file."))
    }
}
comp_info$LFC_cutoff <- as.numeric(comp_info$LFC_cutoff)
if(sum(comp_info$LFC_cutoff<0)>0) stop("'LFC_cutoff' in comparison file is required to be non-negative!")
comp_info$Group_test <- as.character(comp_info$Group_test)
comp_info$Group_ctrl <- as.character(comp_info$Group_ctrl)
if(setDefault){
    A <- cbind(CompareName=rownames(comp_info),comp_info)
    file.rename(config$comparison_file,paste0(config$comparison_file,".bk"))
    write.csv(A,file=config$comparison_file,row.names=F)
    message("-----> comparison file (",basename(config$comparison_file),") is updated with some default values!")
    message("\tThe original comparison file is renamed as ...bk")
}

## read the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,check.names=F,as.is=T)
checkConsistConfigMeta(config,meta)
rownames(meta) <- meta[,config$sample_name]
# use first group name in comparison file for group information
colnames(meta) <- gsub("group","group.org",colnames(meta))
meta <- cbind(group=apply(meta[,config$sample_group,drop=F],1,function(x)return(paste(x,sep="."))),meta)
# set all Group_name to be charactor
for(i in unique(comp_info$Group_name))meta[,i] <- as.character(meta[,i])

## read the gene quantification input ----
message("====== reading gene quantification ...")
estCount <- effeL <- logTPM <- yaxisLab <- NULL
if(!is.null(config$prj_path)){
    estCount <- readData(paste0(config$prj_path,"/combine_rsem_outputs/genes.estcount_table.txt"))
    effeL <- readData(paste0(config$prj_path,"/combine_rsem_outputs/genes.effective_length.txt"))
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
checkSampleName(rownames(meta),colnames(estCount))
if(!is.null(estCount)) estCount <- estCount[,rownames(meta)]
if(!is.null(effeL)) effeL <- effeL[,rownames(meta)]
estCount <- estCount[apply(estCount,1,function(x)return(sum(x>=config$min_count)))>=config$min_sample,]

## covariates removal ----
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
                          prior=config$count_prio)
}
if(is.null(logTPM)){
    stop("Gene quantification is missing, please provide either raw counts with effective length or TPM")
}
logTPM <- logTPM[rownames(estCount),rownames(meta)]
## sample alias ------
if(!is.null(config$sample_alias)){
    rownames(meta) <- colnames(logTPM) <- colnames(estCount) <- meta[,config$sample_alias]
}
saveRDS(estCount,file=paste0(config$output,"/",config$prj_name,"_estCount.rds"))
## gene definition file ---------
message("====== reading gene annotation ...")
if(!is.null(config$gene_annotation)){
    ProteinGeneName <- read.csv(config$gene_annotation,row.names=1,as.is=T)
}else{
    ProteinGeneName <- data.frame(id=rownames(estCount),
                                  UniqueID=rownames(estCount),
                                  Gene.Name=rownames(estCount),
                                  Biotype='unknown')
}
## comparison -----------
message("====== Starting DEG analyses ...")
DEGs <- Batch_DEG(estCount, meta, comp_info,core=config$core)
message("Formating the DEG results")
compRes <- formatQuickOmicsResult(DEGs,logTPM,meta[,"group"],ProteinGeneName)
data_results <- compRes$Dw
results_long <- compRes$Dl
## produce the network for quickOmics ----------
message("====== gene network generation ...")
network <- getNetwork(logTPM,
                      cor_cutoff=config$gene_network_cor_cutoff,
                      p_cutoff=config$gene_network_p_cutoff,
                      variableN=config$gene_network_high_variable_N,
                      edge_max=as.numeric(config$gene_network_max_edge),
                      edge_min=as.numeric(config$gene_network_min_edge),
                      core=config$core)
save(network,file=paste0(config$output,"/",config$prj_name,"_network.RData"))
## save the R object for quickOmics--------
message("====== saving QuickOmics object ...")
data_wide <- logTPM
data_long <- melt(as.matrix(logTPM))
colnames(data_long) <- c("UniqueID","sampleid","expr")
data_long <- cbind(data_long,group=meta[data_long$sampleid,config$sample_group])


MetaData <- formatQuickOmicsMeta(meta,names(DEGs))
save(data_results,results_long,
     data_wide,data_long,
     MetaData,ProteinGeneName,
     comp_info,
     yaxisLab,
     file=paste0(config$output,"/",config$prj_name,".RData"))
## save the project csv file -------
write.csv(data.frame(Name=ifelse(is.null(config$prj_title),
                                 config$prj_name,
                                 config$prj_title),
                     ShortName=config$prj_name,
                     ProjectID=config$prj_name,
                     Species=config$species, 
                     ExpressionUnit=yaxisLab),
          file=paste0(config$output,"/",config$prj_name,".csv"),
          row.names=F,quote=F)

## finishing -----
#system(paste0("cp ",config$output,"/",config$prj_name,"* "))
message("=================================================\nResults are saved in ",config$output)
system(paste0("cp ",config$output,"/",config$prj_name,"* ",sys_config$QuickOmics_path))
#message("-----> to check the results on QuickOmics, execute the following command on ngs server")
#message("\t\tcp ",config$output,"/",config$prj_name,"* /srv/shiny-server/Quickomics/unlisted/.")
#message("\t\t\t And then Please visit: http://ngs.biogen.com:3838/Quickomics/?unlisted=",config$prj_name)
message("\n-----> Please visit: http://ngs.biogen.com:3838/Quickomics/?serverfile=",config$prj_name)

message("\nPowered by the Computational Biology Group [zhengyu.ouyang@biogen.com]")

sink(paste0(config$output,"/session.EArun"))
sessionInfo()
sink()
