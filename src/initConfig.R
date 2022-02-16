
args = commandArgs(trailingOnly=T)
if(length(args)<2){
    stop("A path to a RNAseq project downloaded from DNAnexus is required!")
}
message("loading resource ...")
suppressMessages(source(paste0(args[1],"utility.R"),chdir=T))

configTmp <- yaml::read_yaml(paste0(args[1],"config.tmp.yml"))
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

<<<<<<< HEAD
strMeta <- paste0(strOut,"/sampleMeta.csv")
strMetaFactor <- paste0(strOut,"/sampleMetaFactor.yml")
strGinfo <- paste0(strOut,"/geneAnnotation.csv")
strComp <- paste0(strOut,"/compareInfo.csv")
strAlignQC <- paste0(strOut,"/alignQC.pdf")
strGlength <- paste0(strOut,"/geneLength.pdf")
## format meta information ----
message("Formatting the sample meta information ...")
pInfo <- list()
if(file.exists(paste0(strPath,"/samplesheet.json")))
    pInfo <- rjson::fromJSON(file=paste0(strPath,"/samplesheet.json"))

# Use data.table fread here - problem is if we strip the comment
# line off sample sheet (for reordering for example), this will
# no longer work. comment.char="#" doesn't seem to work for the 
# NGSone file either
sInfo <- data.frame(
			fread(paste0(strPath,"/samplesheet.tsv")),
				  stringsAsFactors=F)
if(is.null(configTmp$sample_name))
    stop("sample_name needs to be defined in in config tmp (contact admin). Default is Sample_Name")
rownames(sInfo) <- sInfo[,configTmp$sample_name]
sInfo <- sInfo[,apply(sInfo,2,function(x)return(sum(!is.na(x))>0))]

qc <- readQC(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"),rownames(sInfo))
qc <- qc[,matchQCnames(qc,config$qc2meta),drop=F]
meta <- merge(sInfo,qc,by="row.names",sort=F)
rownames(meta) <- meta[,1]
meta <- meta[,-1]
## if Concentration & Volume are both present
strkey <- c("Concentration","Volume")
if(sum(strkey%in%colnames(meta))==2){
    meta <- cbind(meta,Vol_Conc=apply(meta[,strkey],1,prod))
}
system(paste("mkdir -p",strOut))
meta <- metaFactor(meta,strMetaFactor)
write.csv(meta,file=strMeta,row.names=F)#
covariates <- c()
for(i in setdiff(colnames(meta),config$notCovariates)){
    if(length(unique(meta[,i]))>1) covariates <- c(covariates,i)
}
## extract effective length ----------
message("Extracting effective length ...")
a <- getEffectLength(strPath)

## gene annotation file ----
message("Create gene annotation ...")
gInfo <- getAnnotation(paste0(strPath,"/config.json"),config$genome_path)
write.csv(gInfo,file=strGinfo)
## alignment QC plots ---------
message("Plot alignment QC ...")
estT <- readData(paste0(strPath,"/combine_rsem_outputs/genes.tpm_table.txt"),
                 rownames(sInfo))
estC <- readData(paste0(strPath,"/combine_rsem_outputs/genes.estcount_table.txt"),
                 rownames(sInfo))
rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
rownames(estC) <- paste(rownames(estC),gInfo[rownames(estC),"Gene.Name"],sep="|")
qc <- readQC(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"),
             rownames(sInfo))
alignQC(estT,qc,strAlignQC,prioQC=config$qc2meta,estC=estC)
## gene length plots -----
if(length(args)>2 && args[3]=="geneLength"){
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
}

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
prjID <- getProjectID(paste0(strPath,"/config.json"))
configTmp <- readLines(paste0(args[1],"config.tmp.yml"))
configTmp <- gsub("initPrjName",prjID,configTmp)
configTmp <- gsub("initPrjTitle",ifelse(is.null(pInfo$Project$Study_Title),
                                        prjID,
                                        pInfo$Project$Study_Title),configTmp)
configTmp <- gsub("initPrjPath",strPath,configTmp)
configTmp <- gsub("initPrjMeta",strMeta,configTmp)
configTmp <- gsub("initPrjFactor",strMetaFactor,configTmp)
configTmp <- gsub("initSpecies",getSpecies(paste0(strPath,"/config.json")),configTmp)
configTmp <- gsub("initGeneAnnotation",strGinfo,configTmp)
configTmp <- gsub("initOutput",strOut,configTmp)
configTmp <- gsub("initCovariates",paste0("[",paste(covariates,collapse=","),"]"),configTmp)
configTmp <- gsub("initPrjComp",strComp,configTmp)

configTmp <- gsub("shinyOne_Title:",paste("shinyOne_Title:",prjID),configTmp)
configTmp <- gsub("shinyOne_Description:",paste("shinyOne_Description:",
                                               ifelse(is.null(pInfo$Project$Study_Title),
                                                      prjID,
                                                      pInfo$Project$Study_Title)),configTmp)

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

message("Powered by the Computational Biology Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]")

=======
checkConfig(configTmp)
pInfo <- checkInputDir(args[2],sysConfig)
pInfo <- appendMeta(pInfo,
                    configTmp$sample_name,
                    sysConfig$qc2meta)
strMsg <- createInit(args[2],readLines(paste0(args[1],"config.tmp.yml")),pInfo)
>>>>>>> publication

finishInit(strMsg)
saveSessionInfo(paste0(strMsg$strOut,"/session.EAinit"),args[1])

