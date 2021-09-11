rm(list=ls())
args = commandArgs(trailingOnly=T)
if(length(args)<1){
    stop("A path to a RNAseq project downloaded from DNAnexus is required!")
}
## initialization -----
message("loading resource ...")
strPath <- normalizePath(args[2])
if(!dir.exists(strPath)){
    stop(paste(strPath, "is not a valid path"))
}
strOut <- paste0(strPath,"/EA",gsub("\\-","",Sys.Date()),"_0")
ix <- 0
while(dir.exists(strOut)){
    ix <- ix+1
    strOut <- paste0(strPath,"/EA",gsub("\\-","",Sys.Date()),"_",ix)
}
if(length(args)>2 & dir.exists(args[3])) strOut <- normalizePath(args[3])

source(paste0(args[1],"gtf.gz_to_gene_info.R"))
source(paste0(args[1],"getAnnotation.R"))
source(paste0(args[1],"extractEffectiveLength.R"))
source(paste0(args[1],"alignQC.R"))
source(paste0(args[1],"readData.R"))
source(paste0(args[1],"lengthQC.R"))
source(paste0(args[1],"metaFactor.R"))
config <- yaml::read_yaml(paste0(args[1],"sys.yml"))
configTmp <- yaml::read_yaml(paste0(args[1],"config.tmp.yml"))

system(paste("mkdir -p",strOut))
strMeta <- paste0(strOut,"/sampleMeta.csv")
strMetaFactor <- paste0(strOut,"/sampleMetaFactor.yml")
strGinfo <- paste0(strOut,"/geneAnnotation.csv")
strComp <- paste0(strOut,"/compareInfo.csv")
strAlignQC <- paste0(strOut,"/alignQC.pdf")
strGlength <- paste0(strOut,"/geneLength.pdf")
## extract effective length ----------
message("Extracting effective length ...")
a <- getEffectLength(strPath)
## format meta information ----
message("Formatting the sample meta information ...")
studyInfo <- tryCatch(
    {
        rjson::fromJSON(gsub('^.','',readLines(paste0(strPath,"/samplesheet.tsv"),n=1)))
    },
    error=function(cond){
        list(Sample_count=NA,#as.numeric(gsub('#Sample number: ','',readLines(paste0(strPath,"/samplesheet.tsv"),n=1))),
             Study_Title=NULL)
    }
)
sInfo <- read.table(paste0(strPath,"/samplesheet.tsv"),sep="\t",
                    header=T,as.is=T,skip=1,comment.char="")
if(is.null(configTmp$sample_name))
    stop("sample_name needs to be defined in in config tmp (contact admin). Default is Sample_Name")
rownames(sInfo) <- sInfo[,configTmp$sample_name]
sInfo <- sInfo[,apply(sInfo,2,function(x)return(sum(!is.na(x))>0))]

qc <- readQC(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"))
qc <- qc[,matchQCnames(qc,config$qc2meta),drop=F]
meta <- merge(sInfo,qc,by="row.names")
rownames(meta) <- meta[,1]
meta <- meta[,-1]
## if Concentration & Volume are both present
strkey <- c("Concentration","Volume")
if(sum(strkey%in%colnames(meta))==2){
    meta <- cbind(meta,Vol_Conc=apply(meta[,strkey],1,prod))
}
meta <- metaFactor(meta,strMetaFactor)
write.csv(meta,file=strMeta,row.names=F)#
covariates <- c()
for(i in setdiff(colnames(meta),config$notCovariates)){
    if(length(unique(meta[,i]))>1) covariates <- c(covariates,i)
}

## gene annotation file ----
message("Create gene annotation ...")
gInfo <- getAnnotation(paste0(strPath,"/config.json"),config$genome_path)
write.csv(gInfo,file=strGinfo)
## alignment QC plots ---------
message("Plot alignment QC ...")
estT <- readData(paste0(strPath,"/combine_rsem_outputs/genes.tpm_table.txt"))
rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
qc <- readQC(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"))
alignQC(estT,qc,strAlignQC,prioQC=config$qc2meta)
## gene length plots -----
message("Plot gene length against expression ...")
estT <- readData(paste0(strPath,"/combine_rsem_outputs/genes.tpm_table.txt"))
estC <- readData(paste0(strPath,"/combine_rsem_outputs/genes.estcount_table.txt"))
effL <- readData(paste0(strPath,"/combine_rsem_outputs/genes.effective_length.txt"))
pdf(strGlength,width=9)
par(mar=c(3,3,2,0.5)+0.1,mgp=c(1.1,0.2,0),tcl=-0.1)
lengthQC(estT,gInfo[rownames(estT),"Length"],main="TPM")
lengthQC(estC,gInfo[rownames(estC),"Length"],main="est Count")
lengthQC(estT,effL[rownames(estT),],main="TPM")
lengthQC(estC,effL[rownames(estC),],main="est Count")
a <- dev.off()

## create an empty comparison file -------
message("Create empty comparison template ...")
comTitle <- c("CompareName",
              "Subsetting_group",
              "Model",
              "Covariate_levels",
              "Group_name",
              "Group_test",
              "Group_ctrl",
              "Analysis_method",
              "Shrink_logFC",
              "LFC_cutoff")
cat(paste(comTitle,collapse=","),"\n",sep="",file=strComp)

## generate a config file ----------
configTmp <- readLines(paste0(args[1],"config.tmp.yml"))
configTmp <- gsub("initPrjName",getProjectID(paste0(strPath,"/config.json")),configTmp)
configTmp <- gsub("initPrjTitle",ifelse(is.null(studyInfo$Study_Title),
                                        getProjectID(paste0(strPath,"/config.json")),
                                        studyInfo$Study_Title),configTmp)
configTmp <- gsub("initPrjPath",strPath,configTmp)
configTmp <- gsub("initPrjMeta",strMeta,configTmp)
configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
configTmp <- gsub("initSpecies",getSpecies(paste0(strPath,"/config.json")),configTmp)
configTmp <- gsub("initGeneAnnotation",strGinfo,configTmp)
configTmp <- gsub("initOutput",strOut,configTmp)
configTmp <- gsub("initCovariates",paste0("[",paste(covariates,collapse=","),"]"),configTmp)
configTmp <- gsub("initPrjComp",strComp,configTmp)

cat(paste(configTmp,collapse="\n"),"\n",sep="",file=paste0(strOut,"/config.yml"))

## finishing --------
message("==========================================")
message("ExpressionAnalysis project folder is created at ",strOut)
message("-----> 'EAqc' can be used to identify the covariates to be adjusted as:")
message("\t\tEAqc ",strOut,"/config.yml\n\n")
message("----->'EArun' can be used to obtain the QuickOmics objects after comparison definition file is updated:")
message("\t\t\t",strComp)
message("\t\tEArun ",strOut,"/config.yml\n\n")
message("-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.\n")

message("Powered by the Computational Biology Group [zhengyu.ouyang@biogen.com]")


sink(paste0(strOut,"/session.EAinit"))
sessionInfo()
sink()

