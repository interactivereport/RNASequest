#args <- c("/home/zouyang/projects/quickOmics/src/","./config.yml")
args = commandArgs(trailingOnly=T)
if(length(args)<2){
    message("'EAinit' can be used to create a config file for an RNAseq project")
    stop("config yaml file is required!")
}
message("loading resource ...")
source(paste0(args[1],"utility.R"),chdir=T)
config <- sapply(yaml::read_yaml(args[2]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

## loading EA data --------
checkConfig(config)
D <- getEAData(config)
D <- useAlias(config,D)

## plotting sequencing QC ----
tmp <- plotAlignQC(2^D$logTPM-config$count_prior,
            paste0(config$output,"/sequenceQC.pdf"),
            estC=D$counts,
            qc=D$seqQC,
            prioQC=sysConfig$qc2meta,
            gInfo=D$gInfo,
            replot=config$seqQC_replot)

## plotting gene length against expression-----
tmp <- plotGeneLength(config,D$counts,effL=D$effLength,logTPM=D$logTPM,gInfo=D$gInfo)

## plot PC analysis ----------
newLogTPM <- plotPCanlaysis(config,D$logTPM,D$meta,estC=D$counts,effL=D$effLength)
if(!is.null(newLogTPM)){
    plotAlignQC(2^newLogTPM-config$count_prior,
                paste0(config$output,"/Adjusted_sequenceQC.pdf"),
                estC=D$counts,
                qc=D$seqQC,
                prioQC=sysConfig$qc2meta,
                gInfo=D$gInfo,
                replot=T)
    
}

## finishing ----
finishQC(list(comparison_file=config$comparison_file,output=config$output))
saveSessionInfo(paste0(config$output,"/session.EAqc"),args[1])
