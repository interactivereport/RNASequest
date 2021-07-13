suppressMessages(require(data.table))
getEffectLength<-function(strPath){
    strF <- paste0(strPath,"/combine_rsem_outputs/genes.effective_length.txt")
    if(file.exists(strF)) return()
    if(dir.exists(paste0(strPath,"/rsem"))){
        D <- NULL
        selColumn <- c("gene_id","effective_length")
        for(i in list.files(paste0(strPath,"/rsem"),"genes.results.gz",full.names=T)){
            X <- fread(i,sep="\t",header=T)[,..selColumn]#,as.is=T,row.names=1
            colnames(X) <- c(selColumn[1],gsub(".genes.results.gz",paste0("|",selColumn[2]),basename(i)))
            if(is.null(D)) D <- X
            else D <- merge(D,X,by=selColumn[1],all=T)
        }
        write.table(D,file=strF,sep="\t",row.names=F,quote=F)
        #return(D)
    }
}