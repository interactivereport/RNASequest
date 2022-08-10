
reSubN <- 2
sbatchDEG <- function(estC,meta,comp_info,StrOutput,srcDir,core){
  compList <- getCompList(comp_info)
  
  strOut <- paste0(StrOutput,"/sbatchOut/")
  system(paste0("rm -rf ",strOut,";mkdir -p ",strOut))
  save(meta,compList,comp_info,srcDir,core,file=paste0(strOut,"env.rdata"))
  
  message("--- Submitting SLURM jobs and monitoring ---")
  sInfo <- list(comNames=names(compList),jID=NULL,runCount=NULL)
  while(length(sInfo$comNames)>0){
    jID <- sInfo$jID
    sInfo <- squeue(strOut,srcDir,sInfo$comNames,core,jID,sInfo$runCount)
    if(is.null(jID)) cat(sInfo$jID)
    Sys.sleep(60)
  }
  DEGs <- list()
  for(one in names(compList)){
    strRDS <- paste0(strOut,one,".rds")
    if(!file.exists(strRDS)){
      message("=====Failed DEG: ",one)
      print(compList[[one]])
      message("==========")
    }
    DEGs <- c(DEGs,readRDS(strRDS))
  }
  return(DEGs)
}
getCompList <- function(comp_info){
  compList <- sapply(group_DEG(comp_info), # function defined in "QuickOmics_DEG.R"
                     function(x)return(list(rownames(x))))
  names(compList) <- paste0("comparisonModel_",1:length(compList))
  return(compList)
}
sbatch <- function(strOut,srcDir,comNames,core,jID=NULL){
  if(is.null(jID)) jID <- paste0("j",sample(10:99,1))
  sapply(comNames,function(one){
    oneScript <- gsub("jID",jID,
                      gsub("jName",one,
                           gsub("wkPath",strOut,
                                gsub("CoreNum",core,
                                     gsub("sysPath",srcDir,
                                          sbatchScript)))))
    cmd <- paste0("Rscript ",srcDir,"/sbatchDEG.R ",strOut," ",one," slurmRun")
    oneScript <- gsub("strCMD",cmd,oneScript)
    strSH <- paste0(strOut,one,".sh")
    cat(oneScript,file=strSH)
    system(paste("sbatch",strSH))
  })
  return(jID)
}
squeue <- function(strOut,srcDir,comNames,core,jID=NULL,runCount=NULL){
  cat(".")
  if(is.null(jID)){
    jID <- sbatch(strOut,srcDir,comNames,core)
    runCount <- setNames(rep(1,length(comNames)),comNames)
    Sys.sleep(1)
  }
  runJobs <- getRunJobs(jID)
  resubJobs <- c()
  for(one in comNames[!comNames%in%runJobs]){
    if(!file.exists(paste0(strOut,one,".rds"))){
      if(runCount[one]<reSubN){
        message("resubmitting ",one)
        resubJobs <- c(resubJobs,one)
      }
      runCount[one] <- runCount[one]+1
    }else{
      system(paste0("cat ",strOut,one,".log"))
    }
  }
  if(length(resubJobs)>0) sbatch(strOut,srcDir,resubJobs,core,jID)
  return(list(comNames=c(runJobs,resubJobs),jID=jID,runCount=runCount))
}
getRunJobs <- function(jID){
  jRec <- system(paste0("squeue | grep '",jID,"_'"),intern=T)
  if(length(jRec)==0) return(NULL)
  
  jNames <- sapply(jRec,function(one){
    oneID <- sapply(strsplit(trimws(gsub("  *"," ",one))," "),head,1)
    tmp <- sapply(strsplit(grep("JobName=",system(paste("scontrol show job",oneID),
                                                  intern=T),
                                value=T),
                           "JobName="),
                  "[[",2)
    return(gsub(paste0("^",jID,"_"),"",sapply(strsplit(tmp," "),head,1)))
  })
  return(jNames)
}
sbatchScript <- "#!/bin/bash
#SBATCH -J jID_jName
#SBATCH -D wkPath
#SBATCH -n CoreNum
#SBATCH -t 72:0:0
#SBATCH -o jName.log
#SBATCH -e jName.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
echo 'end of HOST'
# exit
set -e
env -i bash -c 'source sysPath/.env;eval $condaEnv;strCMD'
echo 'DONE'"

main <- function(args){
  strOut <- args[1]
  comName <- args[2]
  load(paste0(strOut,"env.rdata"))
  estC <- readRDS(list.files(paste0(strOut,".."),"_estCount.rds$",full.names=T)[1])
  source(paste0(srcDir,"/QuickOmics_DEG.R"))
  DEG <- Batch_DEG(estC,meta,comp_info[compList[[comName]],],core=core)
  saveRDS(DEG,file=paste0(strOut,"/",comName,".rds"))
}
args <- commandArgs(trailingOnly=T)
if(length(args)==3){
  print(args)
  message(length(args))
  main(args)
}