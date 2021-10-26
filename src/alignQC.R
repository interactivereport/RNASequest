suppressWarnings(suppressMessages(require(ggplot2)))
suppressWarnings(suppressMessages(require(reshape2)))
alignQC <- function(estT,qc,strPDF,prioQC,topN=c(1,10,30),sIDalias=NULL,estC=NULL){#,50,100
    #estT <- readData(paste0(strPath,"/combine_rsem_outputs/genes.tpm_table.txt"))
    #rownames(estT) <- paste(rownames(estT),gInfo[rownames(estT),"Gene.Name"],sep="|")
    #qc <- readQC(paste0(strPath,"/combine_rnaseqc/combined.metrics.tsv"))

    prioQC <- matchQCnames(qc,prioQC)
    ## if alias provided
    if(!is.null(sIDalias)){
        estT <- estT[,colnames(estT)%in%names(sIDalias)]
        qc <- qc[rownames(qc)%in%names(sIDalias),]
        colnames(estT) <- sIDalias[colnames(estT)]
        rownames(qc) <- sIDalias[rownames(qc)]
    }

    pdfW <- max(nrow(qc)/10+2,6)
    pdf(strPDF,width=pdfW,height=6)
    ## intergenic, intronic and exonic -----
    selN <- c("Exonic_Rate","Intronic_Rate","Intergenic_Rate")
    if(sum(selN%in%colnames(qc))==length(selN)){
        D = melt(as.matrix(qc[,colnames(qc)%in%selN]))
        D$Var1 <- factor(D$Var1,levels=rownames(qc))
        print(ggplot(D,aes(x=Var1,y=value,fill=Var2))+
                  geom_bar(position="stack",stat="identity")+
                  ylab("Fraction of Reads")+xlab("")+
                  ggtitle("Mapped reads allocation")+
                  ylim(0,1)+theme_minimal()+
                  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
                  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                        legend.position = "top")+
                  guides(fill=guide_legend(title="")))
        qc <- qc[,!colnames(qc)%in%selN,drop=F]
    }
    ## top genes ratio----
    topN <- setNames(topN,paste0("Top",topN))
    D <- t(apply(estT,2,function(x){
        x <- sort(x,decreasing=T)
        return(sapply(topN,function(i)return(sum(x[1:i])/sum(x)*100)))
        }))[rownames(qc),,drop=F]
    D <- cbind(sID=factor(rownames(D),levels=rownames(D)),data.frame(D))
    for(i in colnames(D)){
        if(i=="sID") next
        print(ggplot(D,aes_string(x="sID",y=i))+
                  geom_bar(stat="identity")+
                  xlab("")+ylab("% of Total TPM")+
                  ggtitle(paste(i,"genes"))+
                  theme_minimal()+
                  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
                        legend.position = "none"))
    }
    ## union top genes across samples ------
    topUnion <- 30
    p <- plotTopGeneRatio(estT,30)
    write.csv(p$data,file=gsub("pdf","unionTop.TPM.csv",strPDF),row.names=F)
    if(pdfW>8) print(p+theme(aspect.ratio=0.75))
    else print(p)
    # same genes for counts
    if(!is.null(estC)){
        pC <- plotTopGeneRatio(estC,30,selG=levels(p$data[,2]))+
            xlab("% of total counts")+ggtitle(p$labels$title)
        if(pdfW>8) print(pC+theme(aspect.ratio=0.75))
        else print(pC)
    }
    ## all rest qc -----
    qc <- cbind(sID=factor(rownames(qc),levels=rownames(qc)),qc)
    selQC <- colnames(qc)%in%prioQC
    for(i in c(colnames(qc)[selQC],colnames(qc)[!selQC])){
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

plotTopGeneRatio <- function(X,topN,maxN=90,selG=NULL){
    selLevel <- selG
    if(is.null(selG)){
        selG <- unique(as.vector(apply(X,2,function(x)return(names(sort(x,decreasing=T)[1:topN])))))
        while(length(selG)>maxN){
            topN <- topN-1
            selG <- unique(as.vector(apply(X,2,function(x)return(names(sort(x,decreasing=T)[1:topN])))))
        }
    }
    Xsum <- apply(X,2,sum)
    D <- apply(X[selG,],1,function(x)return(100*x/Xsum))
    D <- reshape2::melt(D[,order(apply(D,2,median))])
    if(!is.null(selLevel)) D$Var2 <- factor(D$Var2,selLevel)
    p <- ggplot(D,aes(x=value,y=Var2))+
        geom_point(color="grey50",alpha=0.4,size=1)+
        geom_boxplot(color="#ff7f00",outlier.shape = NA,alpha=0)+
        xlab("% of total TPM")+ylab("")+
        ggtitle(paste("Union of top",topN,"expressed genes"))+
        theme_minimal()+
        theme(axis.text.y=element_text(size=12-(length(selG)-10)/10))
    return(p)
}

#config <- yaml::read_yaml(paste0("/home/zouyang/projects/quickOmics/src/sys.yml"))
#source("/home/zouyang/projects/quickOmics/src/readData.R")
#qc <- readQC(paste0("./combine_rnaseqc/combined.metrics.tsv"))
#matchQCnames(qc,config$qc2meta)


