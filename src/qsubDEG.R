## submit each DEG to clusters

args <- commandArgs(trailingOnly=T)
qsubDEGsrcFile <- "qsubDEGsrc.RData"
qsubDEG <- function(estCount,meta,comp_info,strWK,strSrc,core=8){
    jID <- paste0("j",sample(10:99,1))
    qsubDEGsh <- "#!/bin/bash
#$ -N jID_cName
#$ -wd wkPath
#$ -pe node qsubCore
#$ -l h_rt=50:00:00
#$ -o cName.log
#$ -e cName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH

source /etc/profile.d/modules_bash.sh
module add R/3.5.1
"
    strOut <- paste0(strWK,"/qsubOut/")
    system(paste0("rm -f -R ",strOut,";mkdir -p ",strOut))
    save(estCount,meta,comp_info,core,file=paste0(strOut,qsubDEGsrcFile))
    ## submit each comparison job
    for(i in rownames(comp_info)){
        cat(paste0(gsub("jID",jID,
                        gsub("cName",i,
                             gsub("wkPath",strOut,
                                  gsub("qsubCore",core,qsubDEGsh)))),
                   "Rscript ",strSrc,"qsubDEG.R ",strSrc," ",qsubDEGsrcFile," ",i,"\n"),
            file=paste0(strOut,i,".sh"))
        system(paste0("qsub ",strOut,i,".sh"))
    }
    ## submit monior job depdents on the previous jobs
    message("===== Monitorring all submited DEG jobs ...")
    cat(gsub("wkPath",strOut,
             gsub("qsubCore",1,qsubDEGsh)),
        file=paste0(strOut,"jIDmonitor.sh"))
    system(paste0('qsub -sync y -hold_jid "',jID,'_*" ',strOut,"jIDmonitor.sh"))
    DEGs <- list()
    for(i in rownames(comp_info)){
        strRDS <- paste0(strOut,i,".rds")
        system(paste("cat",gsub("rds$","log",strRDS)))
        if(!file.exists(strRDS)){
            stop(paste("above @",i))
        }
        DEGs <- c(DEGs,readRDS(strRDS))
    }
    return(DEGs)
}

if(length(args)>2 && !grepl("yml$",args[2])){
    suppressMessages(suppressWarnings(source(paste0(args[1],"QuickOmics_DEG.R"))))
    load(args[2])
    DEGs <- Batch_DEG(estCount, meta, comp_info[args[3],],core=core)
    saveRDS(DEGs,file=paste0(args[3],".rds"))
}