
args = commandArgs(trailingOnly=T)
if(length(args)<2){
    message("'EAinit' can be used to create a config file for an RNAseq project")
    stop("config yaml file is required!")
}
message("loading resource ...")
suppressMessages(source(paste0(args[1],"utility.R"),chdir=T))
config <- sapply(yaml::read_yaml(args[2]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

## submit to ShinyOne ------------
Res <- pubShinyOneNEW(c(config,sysConfig))
finishShinyOne(list(shinyApp=sysConfig$shinyApp,ID=Res$id))

##save log
log_file=paste0(args[2], ".EApub.log")
cat(date(), "\n",
    "ShinyOne access: ", sysConfig$shinyApp,"app_project_review.php?ID=",Res$id, "\n\n",
    sep="",  file=log_file )
shinyOneData=Res$shinyOneData
if ("URL" %in% names(shinyOneData)) {
  cat( "URL: ", shinyOneData$URL, "\n", sep="",  file=log_file, append=T )
}
if ("URL_Private" %in% names(shinyOneData)) {
  cat( "Private URL: ", shinyOneData$URL_Private,"\n", sep="",  file=log_file, append=T )
}
cat("\n\nPowered by:",sysConfig$powerby, capture.output(sessionInfo()), sep="\n", file=log_file, append=T)
## R markdown

## finishing
#saveSessionInfo(paste0(config$output,"/session.EApub"),args[1])
