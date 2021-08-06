suppressWarnings(suppressMessages(require(ggplot2)))
suppressWarnings(suppressMessages(require(reshape2)))
alignQC <- function(strPath,gInfo,strPDF,topN=c(1,10,30,50,100)){
    qc <- read.table(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"),
                     sep="\t",header=T,as.is=T,comment.char="",row.names=1,check.names=F,quote="")
    estT <- read.table(paste0(strPath,"/combine_rsem_outputs/genes.tpm_table.txt"),
                       header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    dimnames(qc) <- list(sapply(strsplit(rownames(qc),"_"),function(x)return(gsub(".genome.sorted$","",paste(x[-1],collapse="_")))),
                         gsub("^3","x3",gsub("5","x5",gsub("'","p",gsub(" ","_",gsub("\\%","percentage",colnames(qc)))))))
    qc <- qc[order(rownames(qc)),]
    dimnames(estT) <- list(paste(rownames(estT),gInfo[rownames(estT),'Gene.Name'],sep="|"),
                           sapply(strsplit(sapply(strsplit(colnames(estT),"\\|"),
                                             function(x)return(paste(head(x,-1),
                                                                     collapse="|"))),
                                      "_"),function(x)return(paste(x[-1],
                                                                   collapse="_"))))

    pdf(strPDF,width=max(nrow(qc)/10+2,6))
    ## intergenic, intronic and exonic -----
    selN <- c("Exonic_Rate","Intronic_Rate","Intergenic_Rate")
    D = melt(as.matrix(qc[,colnames(qc)%in%selN]))
    print(ggplot(D,aes(x=Var1,y=value,fill=Var2))+
              geom_bar(position="stack",stat="identity")+
              ylab("Fraction of Reads")+xlab("")+
              ggtitle("Mapped reads allocation")+
              ylim(0,1)+theme_minimal()+
              scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
              theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                    legend.position = "top")+
              guides(fill=guide_legend(title="")))
    qc <- qc[,!colnames(qc)%in%selN]
    ## top genes ratio----
    topN <- setNames(topN,paste0("Top",topN))
    D <- t(apply(estT,2,function(x){
        x <- sort(x,decreasing=T)
        return(sapply(topN,function(i)return(sum(x[1:i])/sum(x)*100)))
        }))[rownames(qc),]
    D <- cbind(sID=rownames(D),data.frame(D))
    for(i in colnames(D)){
        if(i=="sID") next
        print(ggplot(D,aes_string(x="sID",y=i))+
                  geom_bar(stat="identity")+
                  xlab("")+ylab("Percentage")+
                  ggtitle(paste(i,"genes"))+
                  theme_minimal()+
                  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                        legend.position = "none"))
    }
    ## top 30 genes across samples ------
    topG <- unique(as.vector(apply(estT,2,function(x)return(names(sort(x,decreasing=T)[1:30])))))
    estTsum <- apply(estT,2,sum)
    D <- apply(estT[topG,],1,function(x)return(100*x/estTsum))
    D <- melt(D[,order(apply(D,2,median))])
    print(ggplot(D,aes(x=value,y=Var2))+
              geom_boxplot(color="#ff7f00",outlier.shape = NA)+
              geom_point(color="grey50",alpha=0.4,size=1)+
              xlab("% of total TPM")+ylab("")+
              ggtitle("Top expressed genes")+
              theme_minimal())
    ## all rest qc -----
    qc <- cbind(sID=rownames(qc),qc)
    for(i in colnames(qc)){
        if(!is.numeric(qc[1,i])) next
        print(ggplot(qc,aes_string(x="sID",y=i))+
                         geom_bar(stat="identity")+
                         xlab("")+ylab("")+
                         ggtitle(i)+
                         theme_minimal()+
                         theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                               legend.position = "none"))
    }
    ## -----
    a <- dev.off()
}