args = commandArgs(trailingOnly=T)
if(length(args)<2){
    message("'EAinit' can be used to create a config file for an RNAseq project")
    stop("config yaml file is required!")
}
message("loading resource ...")
suppressMessages(source(paste0(args[1],"utility.R"),chdir=T))
config <- sapply(yaml::read_yaml(args[2]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))
config$srcDir <- args[1]
config$ylab <- paste0("log2(TPM+",config$count_prior,")")

## loading EA data --------
checkConfig(config)
D <- getEAData(config)
D <- useAlias(config,D)

## comparison file checking ---------
DEGs <- comparisonAnalysis(config,D$counts,D$meta)

## covariate removal -----
if(!is.null(config$covariates_adjust) && length(config$covariates_adjust)>0){
    batchX <- meta[,config$covariates_adjust,drop=F]
    D$logTPM <- suppressMessages(covariateRM(D$counts,D$effLength,
                                             batchX=batchX,method='limma',
                                             prior=config$count_prior))
    config$ylab <- paste0("log2(adjTPM+",config$count_prior,")")
}

## network analysis -----
saveNetwork(D$logTPM,config)

## formating QuickOmics objects ----
saveQuickOmics(config,D,DEGs)

## finishing ----
finishRun(c(config,sysConfig))
saveSeesionInfo(paste0(config$output,"/session.EArun"),args[1])
