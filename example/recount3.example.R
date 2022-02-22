require(recount3)
#require(biomaRt)

organism <- "mouse"
prjID <- "SRP199678"
strOut <- "./"

## obtain prject through recount3 ---- 
prj <- available_projects(organism=organism)
prj_info <- subset(prj,project==prjID)
rse_gene <- create_rse(prj_info) 

# Extract sample information from sra.sample_attributes ----
meta <- colData(rse_gene)
names(meta) <- gsub("__","_",
                    gsub("[[:punct:]]","_",
                         gsub("%","percent",
                              names(meta))))
grp <- t(sapply(meta$sra_sample_attributes,function(x){
    x <- unlist(strsplit(x,"\\|"))
    return(setNames(trimws(sapply(strsplit(x,";;"),tail,1)),
                    trimws(sapply(strsplit(x,";;"),head,1))))
}))
selIx <- c()
for(one in names(meta)){
    if(length(unique(meta[[one]]))>1) selIx <- c(selIx,one)
    if(is.character(meta[[one]]))
        meta[[one]] <- gsub("[[:punct:]]","_",meta[[one]])
}
meta <- meta[selIx]
meta <- cbind(meta,grp)
grp <- apply(grp,1,paste,collapse="_")
grpN <- sort(table(grp),decreasing = T)
grpN <- setNames(paste0("grp",1:length(grpN)),names(grpN))
meta$group <- as.vector(grpN[grp])

# Extract gene information ----
gInfo <- as.data.frame(rowRanges(rse_gene))
gInfo$id <- 1:nrow(gInfo)
# checked with biomaRt, it is the same as transcript_length, when only one transcript
gLname <- c(grep("bp_length",colnames(gInfo),value=T),grep("score",colnames(gInfo),value=T))[1]
gLength <- matrix(rep(gInfo[,gLname],nrow(meta)),nrow=nrow(gInfo),
                  dimnames=list(rownames(gInfo),rownames(meta)))
colnames(gInfo) <- gsub("gene_id","UniqueID",
                        gsub("gene_name","Gene.Name",
                             gsub("score","Length",
                                  colnames(gInfo))))
# counts -------
X <- transform_counts(rse_gene)

# saving ------
strPath <- paste0(strOut,prjID,"/")
if(!dir.exists(strPath)) dir.create(strPath)

write.csv(meta,file=paste0(strPath,"meta.csv"))
write.csv(gInfo,file=paste0(strPath,"geneAnnotation.csv"))
write.table(gLength,file=paste0(strPath,"geneLength.tsv"),sep="\t",col.names = NA)
write.table(X,file=paste0(strPath,"counts.tsv"),sep="\t",col.names = NA)
cat("metacol,case,ctrl\ngroup,grp2,grp1\n",file=paste0(strPath,"comparison.csv"))
prj_info <- unlist(prj_info)
names(prj_info) <- gsub("n_samples","number_samples",
                        gsub("organism","species",
                             names(prj_info)))
cat(paste(paste(names(prj_info),gsub("[[:punct:]]","_",prj_info),sep=": "),collapse="\n"),"\n",
    sep="",file=paste0(strPath,"project.yml"))



