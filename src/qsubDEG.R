## submit each DEG to clusters

args <- commandArgs(trailingOnly=T)
qsubDEGsrcFile <- "qsubDEGsrc.RData"
qsubDEG <- function(estCount,meta,comp_info,strWK,strSrc,core=8){
    jID <- paste0("j",sample(10:99,1))
    qsubDEGsh <- "#!/bin/bash
#$ -N jID_cName
#$ -wd wkPath
#$ -pe node qsubCore
#$ -l h_rt=1:00:00
#$ -o cName.log
#$ -e cName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH
cat $PE_HOSTFILE
echo 'end of HOST'

source /etc/profile.d/modules_bash.sh
module add R/3.5.1
"
    strOut <- paste0(strWK,"/qsubOut/")
    system(paste0("rm -f -R ",strOut,";mkdir -p ",strOut))
    save(estCount,meta,comp_info,core,file=paste0(strOut,qsubDEGsrcFile))
    ## submit each comparison job
    badNODEs <- qsubSID(rownames(comp_info),jID,strOut,core,qsubDEGsh,strSrc)
    cat(paste(badNODEs,collapse="\n"),file=paste0(strOut,"badNODEs.txt"))
    DEGs <- list()
    for(i in rownames(comp_info)){
        message("===== ",i," =====")
        strRDS <- paste0(strOut,i,".rds")
        system(paste("cat",gsub("rds$","log",strRDS)))
        if(!file.exists(strRDS)){
            message("too many bad nodes: ",strOut,"badNODEs.txt")
            stop(paste("above @",i))
        }
        DEGs <- c(DEGs,readRDS(strRDS))
    }
    return(DEGs)
}
qsubRM <- function(jID,allJOB=F){
    qJob <- tail(system("qstat",intern=T),-2)
    if(!is.null(qJob)){
        qJob<- sapply(strsplit(qJob," "),function(x){
            x <- x[nchar(x)>0]
            if(grepl(jID,x[3])){
                if(x[5]=="Eqw" || allJOB){
                    message("kill job ",x[1],": ",x[3])
                    system(paste("qdel",x[1]))
                    return(F)
                }
                return(T)
            }
            return(F)
        })
    }
    return(sum(qJob))
}
qsubCheckStatus <- function(jID,strOut,qsubDEGsh,sID){
    message("----- Monitorring all submited DEG jobs ...")
    ## check if any job not in good stat
    qsubRM(jID)
    reN <- 0
    maxN <- 60 # maximun 1hr
    while(reN<maxN){
        if(qsubRM(jID)==0)
            break
        Sys.sleep(60) #wait for 1 minutes
        reN <- reN+1
        if(reN>1) message(reN," minutes")
        if(reN>=maxN){
            a <- qsubRM(jID,T)
            break
        }
        # annot use hold_jid since some jobs will enter Eqw state rather than quit
        #strQsub <- paste0(strOut,"jIDmonitor.sh")
        #if(!file.exists(strQsub))
        #    cat(gsub("wkPath",strOut,
        #             gsub("qsubCore",1,
        #                  gsub("h_rt=1:00:00","h_rt=0:05:00",qsubDEGsh))),
        #        file=strQsub)
        #system(paste0('qsub -sync y -hold_jid "',jID,'_*" ',strQsub))
    }

    ix <- sapply(sID,function(x)return(!file.exists(paste0(strOut,x,".rds"))))
    return(ix)
}
qsubSID <- function(sID,jID,strOut,core,qsubDEGsh,strSrc,reN=0,badNodes=NULL){
    for(i in sID){
        strQsub <- paste0(strOut,i,".sh")
        if(!file.exists(strQsub))
            cat(paste0(gsub("jID",jID,
                            gsub("cName",i,
                                 gsub("wkPath",strOut,
                                      gsub("qsubCore",core,qsubDEGsh)))),
                       "Rscript ",strSrc,"qsubDEG.R ",strSrc," ",qsubDEGsrcFile," ",i,"\n"),
                file=strQsub)
        if(is.null(badNodes)){system(paste("qsub",strQsub))}
        else{system(paste0("qsub -l h='!(",paste(badNODEs,collapse="|"),")' ",strQsub))}
    }
    ix <- qsubCheckStatus(jID,strOut,qsubDEGsh,sID)
    if(sum(ix)>0 && reN<5)
        badNodes <- qsubSID(sID[ix],jID,strOut,core,qsubDEGsh,strSrc,reN+1,
                            unique(c(badNodes,getBadNodes(sID[ix],strOut))))
    return(badNodes)
}
getBadNodes <- function(sID,strOut){
    badNODEs <- c()
    for(i in sID){
        strLog <- paste0(strOut,i,".log")
        if(!file.exists(strLog)) next
        tmp <- readLines(strLog)
        badNODEs <- c(badNODEs,tmp[1:(grep("end of HOST",tmp)[1]-1)])
    }
    return(unique(badNODEs))
}

if(length(args)>2 && !grepl("yml$",args[2])){
    suppressMessages(suppressWarnings(source(paste0(args[1],"QuickOmics_DEG.R"))))
    load(args[2])
    DEGs <- Batch_DEG(estCount, meta, comp_info[args[3],],core=core)
    reN <- 0
    while(reN<5 && tryCatch({
        saveRDS(DEGs,file=paste0(args[3],".rds"))
        F
    },error=function(err){
        print(err)
        T
    })){
        Sys.sleep(1)
        reN <- reN+1
    }
}