
args = commandArgs(trailingOnly=T)
if(length(args)<2){
    stop("A path to a RNAseq project downloaded from DNAnexus is required!")
}
message("loading resource ...")
source(paste0(args[1],"utility.R"),chdir=T)
configTmp <- yaml::read_yaml(paste0(args[1],"config.tmp.yml"))
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

checkConfig(configTmp)
pInfo <- checkInputDir(args[2],sysConfig)
pInfo <- appendMeta(pInfo,
                    configTmp$sample_name,
                    sysConfig$qc2meta)
strMsg <- createInit(args[2],readLines(paste0(args[1],"config.tmp.yml")),pInfo)


finishInit(strMsg)
saveSessionInfo(paste0(strMsg$strOut,"/session.EAinit"),args[1])

