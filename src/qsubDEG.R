## submit each DEG to clusters

args <- commandArgs(trailingOnly=T)
qsubDEGsrcFile <- "qsubDEGsrc.RData"
qsubDEG <- function(estCount,meta,comp_info,strWK,strSrc,core=8,qsubTime=180){
    qsubH <- floor(qsubTime/60)
    qsubM <- gsub(" ","0",format(qsubTime%%60,width=2,justify="right"))
    jID <- paste0("j",sample(10:99,1))
    qsubDEGsh <- paste0("#!/bin/bash
#$ -N jID_cName
#$ -wd wkPath
#$ -q short.q
#$ -pe node qsubCore
#$ -l h_rt=",qsubH,":",qsubM,":00
#$ -o cName.log
#$ -e cName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH
cat $PE_HOSTFILE
echo 'end of HOST'

")
    strOut <- paste0(strWK,"/qsubOut/")
    system(paste0("rm -f -R ",strOut,";mkdir -p ",strOut))
    comp_info <- cbind(CompareName=rownames(comp_info),comp_info)
    save(estCount,meta,comp_info,core,file=paste0(strOut,qsubDEGsrcFile))
    ## submit each comparison job
    compList <- getCompList(comp_info)
    badNODEs <- qsubSID(compList,jID,strOut,core,qsubDEGsh,strSrc,qsubTime=qsubTime)
    cat(paste(badNODEs,collapse="\n"),file=paste0(strOut,"badNODEs.txt"))
    DEGs <- list()
    for(i in names(compList)){
        message("===== ",i,": ",paste(compList[i],collapse="; ")," =====")
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
getCompList <- function(comp_info){
	compList <- sapply(group_DEG(comp_info), # function defined in "QuickOmics_DEG.R"
					   function(x)return(list(rownames(x))))
	names(compList) <- paste0("comparisonModel_",1:length(compList))
	return(compList)
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
qsubCheckStatus <- function(jID,strOut,qsubDEGsh,sID,qsubTime){
    message("----- Monitoring all submitted DEG jobs ...")
    ## check if any job not in good stat
    qsubRM(jID)
    reN <- 0
    maxN <- qsubTime # maximun 1hr
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

    ix <- sapply(names(sID),function(x)return(!file.exists(paste0(strOut,x,".rds"))))
    return(ix)
}
qsubSID <- function(sIDs,jID,strOut,core,qsubDEGsh,strSrc,reN=0,badNodes=NULL,qsubTime=180){
    for(i in names(sIDs)){
        Sys.sleep(1)
        strQsub <- paste0(strOut,i,".sh")
        tryIO({
            if(!file.exists(strQsub)){
                cat(paste0(gsub("jID",jID,
                                gsub("cName",i,
                                     gsub("wkPath",strOut,
                                          gsub("qsubCore",core,qsubDEGsh)))),
                           "env -i bash -c 'source ",strSrc,".env;eval $condaEnv;Rscript ",strSrc,"qsubDEG.R ",strSrc," ",qsubDEGsrcFile,' "',paste(sIDs[[i]],collapse=","),'" ',i,"'\n"),
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
    ix <- qsubCheckStatus(jID,strOut,qsubDEGsh,sIDs,qsubTime)
    if(sum(ix)>0 && reN<5)
        badNodes <- qsubSID(sIDs[ix],jID,strOut,core,qsubDEGsh,strSrc,reN+1,
                            unique(c(badNodes,getBadNodes(sIDs[ix],strOut))),
                            qsubTime=qsubTime)
    return(badNodes)
}
getBadNodes <- function(sID,strOut){
    badNODEs <- c()
    for(i in names(sID)){
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
    comNames <- unlist(strsplit(args[3],","))
    DEGs <- Batch_DEG(estCount, meta, comp_info[comNames,],core=core)
    reN <- 0
    tryIO(saveRDS(DEGs,file=paste0(args[4],".rds")))
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
