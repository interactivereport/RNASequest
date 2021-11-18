#args <- c("/home/zouyang/projects/quickOmics/src/","./config.yml")

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
id <- pubShinyOne(config)
finishShinyOne(list(shinyApp=sysConfig$shinyApp,ID=id))

## R markdown

## finishing
saveSeesionInfo(paste0(config$output,"/session.EApub"),args[1])
