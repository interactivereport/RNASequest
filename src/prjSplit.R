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
source(paste0(args[1],"readData.R"))
checkConfig(config)
if(is.null(config$split_meta))
    stop("split_meta is not speicified in the config file")
## read the meta information -----
message("====== reading sample meta information ...")
meta <- read.csv(config$sample_meta,check.names=F,as.is=T)
checkConsistConfigMeta(config,meta)
rownames(meta) <- meta[,config$sample_name]
## read the comparison file ------
message("====== reading comparison information ...")
comp <- read.csv(config$comparison_file,check.names=F,as.is=T)

## create sub projects -------
gInfo <- read.csv(config$gene_annotation,row.names=1,as.is=T)
subProj <- c()
message("Spliting by ",config$split_meta," into: ",paste(unique(meta[,config$split_meta]),collapse=", "))
for(one in unique(meta[,config$split_meta])){
    if(sum(meta[,config$split_meta]==one)<2){
        message("ignore: ",one," with less than 2 samples")
        next
    }
    message("====== Creating sub project: ",one," ...")
    strOut <- paste0(config$output,"/",one)
    subProj <- c(subProj,strOut)
    system(paste("mkdir -p",strOut))
    
    # create the config file
    strConfig <- paste0(strOut,"/config.yml")
    newConfig <- readLines(args[2])
    # replace output and DA output path
    newConfig[grepl("^output",newConfig)] <- paste("output:",strOut)
    newConfig[grepl("^DA_file_outpath",newConfig)] <- paste0("DA_file_outpath: ",strOut,"/DA_Import_Files")

    # create meta file
    strMeta <- paste0(strOut,"/",basename(config$sample_meta))
    oneMeta <- meta[meta[,config$split_meta]==one,]
    write.csv(oneMeta,file=strMeta,row.names=F)
    newConfig[grepl("^sample_meta",newConfig)] <- paste("sample_meta:",strMeta)
    
    # copy comparison definition file
    strCom <- paste0(strOut,"/",basename(config$comparison_file))
    selCom <- gsub(" ","",comp[,"Subsetting_group"])==paste(config$split_meta,one,sep=":")
    if(sum(selCom)>0){
        subCom <- comp[selCom,]
        subCom[,"Subsetting_group"] <- ""
        write.csv(subCom,file=strCom,row.names=F)
    }else{
        write.csv(comp,file=strCom,row.names=F)
    }
    newConfig[grepl("^comparison_file",newConfig)] <- paste("comparison_file:",strCom)

    # change prj_name 
    newConfig[grepl("^prj_name",newConfig)] <- paste0("prj_name: ",config$prj_name,"_",one)
    
    # save new config file
    cat(paste(newConfig,collapse="\n"),"\n",sep="",file=paste0(strOut,"/config.yml"))

    ## the following can be done with EAqc for the subproject by
    ## setting sample_alias to be the same as sample_name
    # plot alignment QC for one sub project
    #message("\tplot QC for the sub-project")
    #alignQC(config$prj_path,
    #        gInfo,
    #        paste0(strOut,"/alignQC.pdf"),
    #        prioQC=sys_config$qc2meta,
    #        sIDalias=setNames(rownames(oneMeta),rownames(oneMeta)))
    
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
message("Powered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")
