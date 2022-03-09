## submit each DEG to clusters

args <- commandArgs(trailingOnly=T)
qsubDEGsrcFile <- "qsubDEGsrc.RData"
qsubDEG <- function(estCount,meta,comp_info,strWK,strSrc,core=8){
    jID <- paste0("j",sample(10:99,1))
    qsubDEGsh <- paste0("#!/bin/bash
#$ -N jID_cName
#$ -wd wkPath
#$ -q short.q
#$ -pe node qsubCore
#$ -l h_rt=1:00:00
#$ -o cName.log
#$ -e cName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH
cat $PE_HOSTFILE
echo 'end of HOST'

",substring(readLines(paste0(strSrc,"sys.yml"),n=1),2),"\n")
    
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
        
        tryIO({system(paste("cat",gsub("rds$","log",strRDS)))})
        one <- NULL
        tryIO({one <- readRDS(strRDS)})
        if(is.null(one)){
            message("too many bad nodes: ",strOut,"badNODEs.txt")
            stop(paste("above @",i))
        }
        DEGs <- c(DEGs,one)
    }
    return(DEGs)
}
qsubRM <- function(jID,allJOB=F){
    qStat <- trySys("qstat",intern=T)
    qJob <- tail(qStat,-2)
    nJob <- length(qJob)
    if(!is.null(qJob) && length(qJob)>0){
        nJob<- sapply(strsplit(qJob," "),function(x){
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
    return(sum(nJob))
}
qsubCheckStatus <- function(jID,strOut,qsubDEGsh,sID){
    message("----- Monitoring all submitted DEG jobs ...")
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
        Sys.sleep(1)
        strQsub <- paste0(strOut,i,".sh")
        tryIO({
            if(!file.exists(strQsub)){
                cat(paste0(gsub("jID",jID,
                                gsub("cName",i,
                                     gsub("wkPath",strOut,
                                          gsub("qsubCore",core,qsubDEGsh)))),
                           "Rscript ",strSrc,"qsubDEG.R ",strSrc," ",qsubDEGsrcFile," ",i,"\n"),
                    file=strQsub)
            }
        })

        if(is.null(badNodes)){
            tryIO({
                a <- system(paste("qsub",strQsub))
                if(a!=0) stop()
                })
        }else{
            tryIO({
                a <- system(paste0("qsub -l h='!(",paste(badNodes,collapse="|"),")' ",strQsub))
                if(a!=0) stop()
                })
        }
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
        tryN <- 0
        tmp <- NULL
        tryIO({tmp <- readLines(strLog)})
        if(is.null(tmp)) next
        badNODEs <- c(badNODEs,sapply(strsplit(tmp[1:(grep("end of HOST",tmp)[1]-1)]," "),head,1))
    }
    return(unique(badNODEs))
}
tryIO <- function(cmd,returnV=NULL,tryN=5){
    N <- 0
    while(N<tryN && tryCatch({
        eval(cmd)
        F
        },error=function(err){
        print(err)
        T
    })){
        Sys.sleep(1)
        N <- N +1
    }
}
trySys <- function(cmd,tryN=5,...){
    N <- 0
    while(N<tryN){
        qStat <- system(cmd,...)
        if(is.null(attr(qStat,"status")) || attr(qStat,"status")==0) break
    }
    return(qStat)
}

if(length(args)>2 && !grepl("yml$",args[2])){
    .libPaths(grep("home",.libPaths(),invert=T,value=T))
    suppressMessages(suppressWarnings(source(paste0(args[1],"QuickOmics_DEG.R"))))
    load(args[2])
    DEGs <- Batch_DEG(estCount, meta, comp_info[args[3],],core=core)
    reN <- 0
    tryIO(saveRDS(DEGs,file=paste0(args[3],".rds")))
    #while(reN<5 && tryCatch({
    #    saveRDS(DEGs,file=paste0(args[3],".rds"))
    #    F
    #},error=function(err){
    #    print(err)
    #    T
    #})){
    #    Sys.sleep(1)
    #    reN <- reN+1
    #}
}
