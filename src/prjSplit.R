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
source(paste0(args[1],"infoCheck.R"))
source(paste0(args[1],"alignQC.R"))
checkConfig(config)
if(is.null(config$split_meta))
    stop("split_meta is not speicified in the config file")
## read the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,check.names=F,as.is=T)
checkConsistConfigMeta(config,meta)
rownames(meta) <- meta[,config$sample_name]
## create sub projects -------
gInfo <- read.csv(config$gene_annotation,row.names=1,as.is=T)
subProj <- c()
for(one in unique(meta[,config$split_meta])){
    if(sum(meta[,config$split_meta]==one)<2){
        message("ignore: ",one," with less than 2 samples")
        next
    }
    message("====== Creating sub project: ",one," ...")
    strOut <- paste0(config$output,"/one")
    subProj <- c(subProj,strOut)
    system(paste("mkdir -p",strOut))
    
    # create the config file
    strConfig <- paste0(strOut,"/config.yml")
    system(paste("cp",args[2],strConfig))
    system(paste0('sed -i "s|^output.*|output: ',strOut,'|" ',strConfig))
    
    # create meta file
    strMeta <- paste0(strOut,"sampleMeta.csv")
    write.csv(meta[meta[,config$split_meta]==one,],file=strMeta)
    system(paste0('sed -i "s|^sample_meta.*|output: ',strMeta,'|" ',strConfig))

    # copy comparison definition file
    strCom <- paste0(strOut,"/compareInfo.csv")
    system(paste("cp",config$comparison_file,strCom))
    system(paste0('sed -i "s|^comparison_file.*|output: ',strCom,'|" ',strConfig))
    
    # plot alignment QC for one sub project
    alignQC(config$prj_path,
            gInfo,
            paste0(strOut,"/alignQC.pdf"),
            prioQC=sys_config$qc2meta,
            sIDalias=setNames(rownames(meta),rownames(meta)))
    
    sink(paste0(strOut,"/session.EAsplit"))
    sessionInfo()
    sink()
    
}
## finishing -----
message("==========================================")
message("Several ExpressionAnalysis project folders are created at:")
message(paste(subProj,collapse="\n"))
message("-----> 'EAqc' can be used to identify the covariates for each sub-project by using config in the folder:")
message("-----> 'EArun' can be used to obtain the QuickOmics objects for each sub-project after comparison definition file is updated:")
message("")
message("Powered by the Computational Biology Group [zhengyu.ouyang@biogen.com]")
