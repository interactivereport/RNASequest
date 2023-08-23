args = commandArgs(trailingOnly=T)
if(length(args)<2){
    #message("'EArun can be used to create a config file for an RNAseq project")
    stop("config yaml file is required!")
}
message("loading resource ...")
source(paste0(args[1],"utility.R"),chdir=T)
config <- sapply(yaml::read_yaml(args[2]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))
config$srcDir <- args[1]
config$ylab <- paste0("log2(TPM+",config$count_prior,")")

## loading EA data --------
config <- checkConfig(config)
a <- checkShinyTestSetting(sysConfig)
D <- getEAData(config,withCom=T)
D <- useAlias(config,D)
D <- lowCountFiltering(config,D)

## check correlation between covariates and comparison groups -----
a <- plotCovBio(config,D$meta,D$comp_info)

## comparison analysis ---------
DEGs <- comparisonAnalysis(config,D$counts,D$meta,D$comp_info)
plotDEG_MA(DEGs,D$logTPM,D$meta,D$comp_info,config)

## covariate removal -----
if(!is.null(config$covariates_adjust) && length(config$covariates_adjust)>0){
    batchX <- D$meta[,config$covariates_adjust,drop=F]
    D$logTPM <- suppressMessages(covariateRM(D$counts,D$effLength,
                                             batchX=batchX,method=config$covariates_method,#'limma',
                                             prior=config$count_prior))
    config$ylab <- paste0("log2(adjTPM+",config$count_prior,")")
}

## network analysis -----
saveNetwork(D$logTPM,config)

## formating QuickOmics objects ----
saveQuickOmics(config,D,DEGs)

## finishing ----
finishRun(c(config,sysConfig))
saveSessionInfo(paste0(config$output,"/session.EArun"),args[1])
