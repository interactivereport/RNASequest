
lengthQC <- function(X,gL,...){
    if(!is.null(ncol(gL)))
        plotGeneLength(cbind(Length=apply(gL,1,median),X),list(mean=mean,median=median,SD=sd),
                       xlab="log2(median gene effective length +1)",...)
    else
        plotGeneLength(cbind(Length=gL,X),list(mean=mean,median=median,SD=sd),
                       xlab="log2(GTF gene length +1)",...)
}

plotGeneLength <- function(X,funs,...){
    imageCOL <- c("#FFFFFFFF","#3300FF","#2D1CFF","#2839FF","#2255FF","#1C71FF","#178EFF","#11AAFF",
                  "#0BC6FF","#06E3FF","#00FFFF","#00FFFF","#17FFE3","#2DFFC6","#44FFAA","#5BFF8E",
                  "#71FF71","#88FF55","#9FFF39","#B5FF1C","#CCFF00",
                  "#CCFF00","#D2F400","#D7E800","#DDDD00","#E3D200","#E8C600","#EEBB00","#F4B000",
                  "#F9A400","#FF9900","#FF9900","#F68800","#EC7700","#E36600","#D95500","#D04400",
                  "#C63300","#BD2200","#B31100","#AA0000")
    X <- log2(X+1)
    for(i in names(funs)){
        y <- apply(X[,-1],1,funs[[i]])
        x <- X[,"Length"]
        #x <- X[y>0,"Length"]
        #y <- y[y>0]
        plot(c(),c(),xlim=range(x),ylim=range(y),
             ylab=paste(i,"of log2(Expression+1)"),...)
        index <- rep(T,length(x))
        tryM <- try(f1 <- MASS::kde2d(x,y,n=c(90,60)),silent=T)
        if(!is.null(names(tryM))){
            image(f1,col=imageCOL,add=T)
            imageZero <- diff(range(f1$z))/(length(imageCOL)-1) # ~2% dots
            index <- apply(cbind(x,y),1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
        }
        points(x[index],y[index],pch=20,col=imageCOL[2],cex=1)
    }
}