rm(list=ls())
.libPaths(grep("home",.libPaths(),invert=T,value=T))
args = commandArgs(trailingOnly=T)
if(length(args)<1){
    stop("A path to a RNAseq project downloaded from DNAnexus is required!")
}

config <- sapply(yaml::read_yaml(args[2]),unlist)
sys_config <- yaml::read_yaml(paste0(args[1],"sys.yml"))

## submit to ShinyOne ------------
strF <- paste0(sys_config$QuickOmics_publish_folder,config$prj_name,".RData")
if(file.exists(strF)){
    stop("The project already exists in ShinyOne!\nPlease remove the record and associated files or change prj_name and re-EArun!")
}

shinyOneData <- config[grep("^shinyOne_",names(config))]
names(shinyOneData) <- gsub("shinyOne_","",names(shinyOneData))
shinyOneData[['Title']] <- paste(config$prj_name,config$prj_prj_title,sep=": ")
shinyOneData[["Species"]] <- config$species
shinyOneData[["TSTID"]] <- ifelse(grepl("^TST",config$prj_name),substr(config$prj_name,1,8),"")
shinyOneData[["Data_Generated_By"]] <- system("whoami",intern=T)
shinyOneData[["Date"]] <- as.character(Sys.Date())
shinyOneData[["URL"]] <- paste0(sys_config$QuiclOmics_publish_link,config$prj_name)

shinyOneCMD <- paste0("curl -s -k -X POST -d 'data={",
                      paste(paste0('"',names(shinyOneData),'": "',shinyOneData,'"'),collapse = ", "),
                      "}' '",sys_config$shinyApp,"api_add_project.php?api_key=lnpJMJ5ClbuHCylWqfBY8BoxxdrpU0'")
# for real job
res <- system(shinyOneCMD,intern=T)
# debug
#res <- "{\"Status\":true,\"ID\":148}"

shinyMsg <- tryCatch({
    rjson::fromJSON(res)
},error=function(eMsg){
    stop(paste0(paste(res,collapse="\n"),
                "\nPlease contact Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]"))
})
if(!shinyMsg$Status){
    stop(paste0(paste(paste(names(shinyMsg),shinyMsg,sep=":"),collapse="\n"),
                "\nPlease contact Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]"))
}
system(paste0("cp ",config$output,"/",config$prj_name,"* ",sys_config$QuickOmics_publish_folder))

message("=================================================\nShinyOne access: ",
        sys_config$shinyApp,"app_project_review.php?ID=",shinyMsg$ID)
message("\nPowered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")


## create publish folders --------
strOut <- paste0(config$output,"/publish/")
strData <- paste0(strOut,"data/")
system(paste("mkdir -p",strData))

## copy raw DNAnexus files (counts and effective length)----
system(paste0("mkdir -p ",strData,"combine_rsem_outputs"))
system(paste0("cp ",config$prj_path,"/combine_rsem_outputs/genes.estcount_table.txt ",strData,"combine_rsem_outputs/."))
system(paste0("cp ",config$prj_path,"/combine_rsem_outputs/genes.effective_length.txt ",strData,"combine_rsem_outputs/."))

## readin process script -----
runScript <- readLines(paste0(args[1],"quickOmics.R"))

## create utility functions -----
system(paste0("echo '' > ",strOut,"utility.R"))
for(one in sapply(strsplit(grep("^source",runScript,value=T),"\"|\'"),function(x)return(grep(".R$",x,value=T)))){
    if(grepl("qsub",one)) next
    system(paste0("cat ",args[1],one," >> ",strOut,"utility.R"))
}

## write the process script for EApub sections----
secIndex <- grep("## .*----$",runScript)
secSections <- data.frame(start=secIndex,end=c(secIndex[-1]-1,length(runScript)))

pubScript <- unlist(apply(secSections,1,function(x){
    oneScript <- ""
    if(grepl("^## EApub",runScript[x[1]])){
        oneScript <- runScript[x[1]:x[2]]
        if(grepl("^## EApub if",oneScript[1])){
            ind <- grep("^if|^}",oneScript)
            if(length(ind)!=3) stop(paste0("Section error for",oneScript[1]))
            oneScript <- oneScript[-c(ind[1],ind[2]:ind[3])]
        }else if(grepl("^## EApub else",oneScript[1])){
            ind <- grep("^if|^}",oneScript)
            if(length(ind)!=3) stop(paste0("Section error for",oneScript[1]))
            oneScript <- oneScript[-c(ind[1]:ind[2],ind[3])]
        }
    }
    return(oneScript)
}))

cat(paste(c('source("utility.R")',
            'config <- sapply(yaml::read_yaml("data/config.yml"),unlist)',
            pubScript),collapse="\n"),
    file=paste0(strOut,"process.R"))

## prepare/save the configer used files -----
usedConfig <- unique(unlist(sapply(pubScript,function(oneScript){
    return(trimws(sapply(strsplit(strsplit(oneScript,"config\\$")[[1]][-1],",|)|]| |>|<"),head,1)))
})))
newConfig <- config[usedConfig]
names(newConfig) <- usedConfig

for(one in usedConfig){
    if(one=="prj_path"){
        newConfig[[one]] <- "data/"
    }else if(one=="output"){
        newConfig[[one]] <- "./"
    }else if(is.character(newConfig[[one]]) && file.exists(newConfig[[one]])){
        if(!file.exists(paste0(strData,basename(newConfig[[one]])))){
            system(paste("cp",newConfig[[one]],strData))
        }
        newConfig[[one]] <- paste0('data/',basename(newConfig[[one]]))
    }
}

cat(paste(paste0(names(newConfig),": ",
                 sapply(newConfig,function(x){
                     if(length(x)==0) return("")
                     if(length(x)==1)return(x)
                     return(paste0('["',paste(x,collapse='","'),'"]'))
                 })),
          collapse="\n"),
    "\n",sep="",file=paste0(strData,"config.yml"))

## -----
message("The publishable processing code and data is saved in:\n\t",strOut)
message("\nPowered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")

