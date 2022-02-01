suppressMessages(require(recount))

strGSE <- "GSE32465"
strOut <- "./"

strPath <- paste0(strOut,strGSE,"/")
dir.create(strPath)

project_info <- abstract_search(strGSE)
download_study(project_info$project)
load(file.path(project_info$project, "rse_gene.Rdata"))

# Extract sample information from characteristics
meta <- colData(rse_gene)
grp <- sapply(meta$characteristics,function(x)
    return(setNames(trimws(sapply(strsplit(x,":"),tail,1)),
                    trimws(sapply(strsplit(x,":"),head,1))))
    )
meta <- cbind(meta[,-ncol(meta)],t(grp))
grp <- apply(grp,2,paste,collapse="_")
grpN <- sort(table(grp),decreasing = T)
grpN <- setNames(paste0("grp",1:length(grpN)),names(grpN))
meta$group <- grpN[grp]

# Extract gene information
gInfo <- rowData(rse_gene)
gLength <- gInfo[,grepl("length",colnames(gInfo))]
gLength <- matrix(rep(gLength,nrow(meta)),nrow=length(gLength),
                  dimnames=list(names(gLength),rownames(meta)))
gInfo <- cbind(gInfo[,!colnames(gInfo)%in%c("gene_id","symbol"),drop=F],
               id=1:nrow(gInfo),
               UniqueID=gInfo$gene_id,
               Gene.Name=sapply(gInfo$symbol,head,1))
# counts
X <- transform_counts(rse_gene)

# saving 
write.csv(meta,file=paste0(strPath,"meta.csv"))
write.csv(gInfo,file=paste0(strPath,"geneAnnotation.csv"))
write.table(gLength,file=paste0(strPath,"geneLength.tsv"),sep="\t",col.names = NA)
write.table(X,file=paste0(strPath,"counts.tsv"),sep="\t",col.names = NA)
cat("metacol,case,ctrl\ngroup,grp2,grp1\n",file=paste0(strPath,"comparison.csv"))
cat(paste(paste(names(project_info),gsub("[[:punct:]]","_",project_info),sep=": "),collapse="\n"),"\n",
    sep="",file=paste0(strPath,"project.yml"))










