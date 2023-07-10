#standalone functions to create plots for R markdown use

options(stringsAsFactors=F)
options(ggrepel.max.overlaps = Inf)

suppressPackageStartupMessages({
  library(plotly)
  library(reshape2)
  library(tidyverse)
  library(gplots)
  library(ggpubr)
  library(gridExtra)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
  library(car)
  library(colourpicker)
  library(VennDiagram)
  library(factoextra)
  library(Mfuzz)
  library(openxlsx)
  library(visNetwork)
  library(cowplot)
  library(circlize)
  library(ComplexHeatmap)
  library(svglite)
  library(pathview)
  library(ggrastr)
  library(biomaRt)
  library(networkD3)
  library(Hmisc)
  library(stringi)
  library(plyr)
  library(ggExtra)
  library(png)
  library(psych)
  library(broom)
})


sub_samples <- function(data_wide, MetaData, input_samples=NULL) {
  #Only use selected samples
  if (!is.null(input_samples)){
    sel_sample_order=match(input_samples, MetaData$sampleid)
    sel_sample_order=sel_sample_order[!is.na(sel_sample_order)]
    MetaData=MetaData[sel_sample_order, ] #user input sample order
  }
  tmp_sampleid = MetaData$sampleid #reorder data_Wide, also select subset
  data_wide  <- data_wide[apply(data_wide, 1, function(x) sum(length(which(x==0 | is.na(x)))) < 3),]
  sel_sample_order2=match(tmp_sampleid, colnames(data_wide))
  sel_sample_order2=sel_sample_order2[!is.na(sel_sample_order2)]
  data_wide = data_wide[, sel_sample_order2] %>% as.matrix()
  return(list(MetaData=MetaData, data_wide=data_wide))
}

##Make PCA plot
PCA_Plot <- function(data_wide, MetaData, input_samples=NULL, pcnum=c(1, 2), label_samples=TRUE,
                      PCAcolorby="group",  PCAshapeby="none", PCAsizeby="none", PCA_label="sampleid",
                      ellipsoid=FALSE, mean_point = FALSE, rug = FALSE,  PCAcolpalette="Dark2",
                      PCAdotsize=4, PCAfontsize=10 ) {
  
  sub_data=sub_samples(data_wide, MetaData, input_samples)
  MetaData=sub_data[["MetaData"]]
  tmp_sampleid = MetaData$sampleid
  tmp_data_wide=sub_data[["data_wide"]]
 
  tmp_data_wide[is.na(tmp_data_wide)] <- 0 
  pca <- 	prcomp(t(tmp_data_wide),rank. = 10, scale = FALSE)
  percentVar <- 	round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
  scores <- as.data.frame(pca$x)
  rownames(scores) <- tmp_sampleid
  attributes=setdiff(colnames(MetaData), c("Order", "ComparePairs") )
  colsel=match(attributes, colnames(MetaData) )
  scores=cbind(scores, MetaData[, colsel, drop=F])
  
  samples=scores$sampleid
  
  xlabel <- paste("PC",pcnum[1],"(",round(percentVar[pcnum[1]]),"%)",sep="")
  ylabel <- paste("PC",pcnum[2],"(",round(percentVar[pcnum[2]]),"%)",sep="")
  
  PC1 <- paste("PC",pcnum[1],sep="")
  PC2 <- paste("PC",pcnum[2],sep="")
  
  n <- length(unique(as.character(unlist(scores[, colnames(scores)==PCAcolorby]))))
  colorpal = colorRampPalette(brewer.pal(8, PCAcolpalette))(n)
  
  
  if (label_samples==FALSE ) {labels=NULL
  } else {
    label_sel=match(PCA_label, names(scores))
    labels=unlist(scores[, label_sel])	
  }
  
  if (PCAshapeby=="none") {shape_by=19} else {shape_by=PCAshapeby}
  if (PCAsizeby=="none") {size_by=PCAdotsize} else {size_by=PCAsizeby}	
  if (is.numeric(scores[[PCAcolorby]])) {  #when colorby is numeric, don't use color palette
    p <- ggpubr::ggscatter(scores,x =PC1, y=PC2, color =PCAcolorby, shape=shape_by, size =size_by , ellipse = ellipsoid, mean.point = mean_point, rug = rug,
                           label =labels, font.label = PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
  } else {
    p <- ggpubr::ggscatter(scores,x =PC1, y=PC2, color =PCAcolorby, shape=shape_by, size =size_by , palette= colorpal, ellipse = ellipsoid, mean.point = mean_point, rug = rug,
                           label =labels, font.label = PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
  }
  
  p <- ggpubr::ggpar(p, xlab = xlabel, ylab = ylabel)
   p <- p + guides(color = guide_legend(override.aes = list(label="")))
  return(p)
}


#Compute Covariates
source("PC_Covariates.R")
PC_covariates_out <-  function(data_wide, MetaData, input_samples=NULL,select_covariates=NULL,  PC_cutoff=5, FDR_cutoff=0.1, N_col=2){
  sub_data=sub_samples(data_wide, MetaData, input_samples)
  MetaData=sub_data[["MetaData"]]
  tmp_sampleid = MetaData$sampleid
  tmp_data_wide=sub_data[["data_wide"]]
  meta=MetaData[, !(colnames(MetaData) %in% c("sampleid", "Order", "ComparePairs")), drop=FALSE]
  if (!is.null(select_covariates)) {
    meta=meta[, (colnames(meta) %in% select_covariates), drop=FALSE]
  }
  rownames(meta)=MetaData$sampleid
  #browser() #debug
  res<-Covariate_PC_Analysis(tmp_data_wide, meta, out_prefix=NULL, PC_cutoff=PC_cutoff, FDR_cutoff=FDR_cutoff, N_col=N_col)
  return(res)
}


PCA_sig_covariates<-function(res) {
  data.all=res$data.all
  selVar_All=res$selVar_All
  PC_info<-res$PC_info
  selVar_All<-selVar_All%>%arrange(FDR) #sort by FDR
  PCA_plots=vector(mode="list", length=nrow(selVar_All))
  for (i in 1:nrow(selVar_All)) {
        x0=selVar_All$PC[i]
        if (x0=="PC1") {y0="PC2"} else {y0=x0; x0="PC1"}
        x=sym(x0); y=sym(y0); color_by=sym(selVar_All$Covariate[i])
        p<-ggplot(data.all, aes(x=!!x, y=!!y, col=!!color_by))+geom_point()+
          labs(x=PC_info$PC_new[PC_info$PC==x0], y=PC_info$PC_new[PC_info$PC==y0])+ theme_half_open()
        #additional text to add  
        p<-add_sub(p, selVar_All$Significance[i], x=0.2, hjust=0)
        PCA_plots[[i]]=p
        names(PCA_plots)[i]=str_c(selVar_All$PC[i], " vs. ", selVar_All$Covariate[i])
  }
  return(PCA_plots)
}

#plot one covariate vs. one PC. Will choose scatter plot or boxplot based on the covariate type
covariate_vs_PC_plot<-function(res, pc, var, add_text=TRUE) {
  data.all=res$data.all
  selVar=res$selVar_All
  if (!(var %in% names(data.all))) {cat("Covariate ", var, " not in MetaData. Please check the spelling of covariate.\n", sep=""); return(NULL)}
  if (!(pc %in% names(data.all))) {cat(pc, " not in principle component scores. Please check the spelling.\n", sep=""); return(NULL)}
  Num_names=names(select_if(data.all, is.numeric))
  if (var %in% Num_names) { #numeric covariate
    p<-ggplot(data.all, aes(x=!!sym(var), y=!!sym(pc)) )+geom_point()+
      stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+theme_half_open()
  } else {
    p<-ggplot(data.all,  aes(x=!!sym(var), y=!!sym(pc)) )+geom_boxplot()+geom_jitter(alpha=0.7, width=0.1)+theme_half_open() +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  }
  if (add_text) {
    info<-selVar%>%dplyr::filter(PC==pc, Covariate==var)
    if ( nrow(info)>0) {text_info=info$Significance[1]} else {text_info=str_c(var, " vs. ", pc, " not significant.")}
    p<-add_sub(p, text_info, x=0.2, hjust=0)
  }
  return(p)
}



Volcano_Plot <- function(results_long, ProteinGeneName, comparison, volcano_FCcut=1.2, volcano_pvalcut=0.05, volcano_psel = "Padj", 
                         volcano_genelabel="Gene.Name", label_DEG=TRUE, Max_Pvalue=0, Max_logFC=0, Ngenes=50, rasterize="Yes",
                         vlegendpos="bottom", lfontsize=4, yfontsize=14, volcano_label=TRUE) {
  results_long <-
    results_long %>% mutate_if(is.factor, as.character)  %>% left_join(ProteinGeneName, by = "UniqueID")
  test_sel = comparison
  FCcut = log2(as.numeric(volcano_FCcut))
  FCcut_rd=round(FCcut*1000)/1000
  pvalcut = as.numeric(volcano_pvalcut)
  res = results_long %>% filter(test==test_sel) %>%
    filter(!is.na(P.Value)) %>%
    dplyr::mutate (color="Not Significant") %>% as.data.frame() 
  
  
  res$labelgeneid = res[,match(volcano_genelabel,colnames(res))]
  
  if (volcano_psel == "Padj") {
    res$color[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = paste0("Padj","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$Adj.P.Value<pvalcut)==1)] =  paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("Padj","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (Max_Pvalue>0) {
      res<-res%>%mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-Max_Pvalue) ))
    }
  } else { 
    res$color[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = paste0("pval","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$P.Value<pvalcut)==1)] =  paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("pval","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (Max_Pvalue>0) {
      res<-res%>%mutate(P.Value=pmax(P.Value, 10^(0-Max_Pvalue) ))
    }
  }
  res$logFC_ori=res$logFC
  if (Max_logFC>0) {
    res<-res%>%mutate(logFC=ifelse(logFC>=0, pmin(Max_logFC, logFC), pmax(0-Max_logFC, logFC) ) )
  }
 
  if (volcano_psel == "Padj") {
    p <- ggplot(res, aes(x = logFC, y = -log10(Adj.P.Value)))
    ylab <- "-log10(Padj.Value)"
    
    filterSig <- paste0("Padj", "<", pvalcut, " & abs(log2FC)>", FCcut_rd)
    data.label <- filter(res, color == filterSig)
    if (nrow(data.label) > Ngenes) {
      data.label <- top_n(data.label, Ngenes, abs(logFC_ori))
    }
    
  } else {
    filterSig <- paste0("pval", "<", pvalcut, " & abs(log2FC)>",FCcut_rd)
    data.label <- filter(res, color == filterSig)
    if (nrow(data.label) > Ngenes) {
      data.label <- top_n(data.label, Ngenes, abs(logFC_ori))
    }
    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value)))
    ylab <- "-log10(P.Value)"
  }
  
  
  p <- p	+
    scale_color_manual(values = c("grey", "green2","red2"))
  if (rasterize=="Yes") { p<-p+geom_point_rast(aes(color = color), size=0.7, alpha=0.6, na.rm=TRUE, dev="ragg") 
  } else {p<-p+geom_point(aes(color = color), size=0.7) }
  p <- p+
    theme_bw(base_size = 20) +
    geom_hline(yintercept = -log10(pvalcut), colour="grey") +
    geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
    ylab(ylab) + xlab("log2 Fold Change") +
    ggtitle(test_sel) +
    theme(legend.position = vlegendpos, legend.text=element_text(size=yfontsize))
  if (volcano_label) {
    p=p+geom_text_repel(data = data.label,  aes(label=labelgeneid),	size = lfontsize,	box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  }
  p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))
  return(p)
  
}


##Heatmap
#heatmap_subset, #options, random, variable, subset (DEG), gene_list
#maxgenes=1000, heatmap_N_genes (number of genes to lable)
#if gene_list is not NULL, ignore heatmap_subset
#subset, must have comparison name, 
#dendrogram", "Apply Clustering:", c("both" ,"none", "row", "column")
Heatmap_Plot<-function(data_wide, MetaData, results_long, ProteinGeneName, 
                       input_samples=NULL, heatmap_subset="variable", gene_list=NULL, comparison=NULL, maxgenes=1000,
                       heatmap_fccut=1.2, heatmap_pvalcut=0.05, heatmap_psel = "Padj",heatmap_annot="group",dendrogram="both",
                       scale="row", hxfontsizep=10, hyfontsizep=7, heatmap_label="Gene.Name", heatmap_N_genes=100, 
                       heatmap_row_dend=TRUE, heatmap_col_dend=TRUE, cutreerows=0, cutreecols=0, exp_unit="log2(TPM+1)",
                       lowColor="blue", midColor="white", highColor="red", distanceMethod="euclidean", agglomerationMethod="complete",
                       heatmap_row_title="", heatmap_row_title_font_size=16,
                       heatmap_column_title="", heatmap_column_title_font_size=16) {
  sub_data=sub_samples(data_wide, MetaData, input_samples)
  MetaData=sub_data[["MetaData"]]
  tmp_sampleid = MetaData$sampleid
  tmpdat=sub_data[["data_wide"]]
  tmpdat[is.na(tmpdat)] <- 0
  annotation = MetaData
  heatmap_fccut =log2(heatmap_fccut)
  if (heatmap_subset == "subset") {
    #from DEG
    if (is.null(comparison)) {stop("Please provide a comparions name.")} else {heatmap_test=comparison} 
    if (heatmap_psel == "Padj") {
      filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & Adj.P.Value < heatmap_pvalcut) %>%
        dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    } else {
      filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & P.Value < heatmap_pvalcut) %>%
        dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    }
    if(length(filteredGene)>2) {
      gene_list=intersect(rownames(tmpdat), filteredGene ) #user only genes that are in expression.
      tmpdat  <-  tmpdat[gene_list,]
    } else {stop("Not enough genes pass DEG cutoff!")}
  } else if (heatmap_subset == "random") {
      tmpdat=tmpdat[sample(1:nrow(tmpdat), maxgenes),] #this will keep rownames
        #tmpdat <- tmpdat %>% sample_n(maxgenes) #this will remove rownames
  } else if  (heatmap_subset == "variable") {
        dataSD=apply(tmpdat, 1, function(x) sd(x,na.rm=T))
        dataM=rowMeans(tmpdat)
        diff=dataSD/(dataM+median(dataM)) #SD/mean, added median value to penalized lowly expressed genes
        tmpdat=tmpdat[order(diff, decreasing=TRUE)[1:maxgenes], ]	    
  } else if (heatmap_subset == "gene_list") {
    if (is.null(gene_list)) {stop("Please enter a gene list to plot heatmap, or choose another heatmap_subset option.")}
    heatmap_list <- unique(gene_list)
    uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% heatmap_list) | (Protein.ID %in% heatmap_list) | (toupper(Gene.Name) %in% toupper(heatmap_list)))  %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()

    #restore order of the input list
    sel1=match(uploadlist, ProteinGeneName$UniqueID)
    ID_order<-ProteinGeneName[sel1, ]%>%mutate(N1=match(UniqueID, heatmap_list), N2=match(Protein.ID, heatmap_list), 
                                               N3=match(toupper(Gene.Name), toupper(heatmap_list)), N=pmin(N1, N2, N3, na.rm=T))%>%arrange(N)
    tmpdat  <-  tmpdat[ID_order$UniqueID,]
  }
  
  if (nrow(tmpdat)>5000 ) {tmpdat=tmpdat[sample(1:nrow(tmpdat), 5000),]; cat("Reduce data pionts to 5K\n")} #Use at most 5000 genes so the App won't crash
  
  
  df <- data.matrix(tmpdat)
  #use selected gene label
  sel=match(rownames(df), ProteinGeneName$UniqueID)
  selCol=match(heatmap_label, names(ProteinGeneName))
  
  if (sum(is.na(sel))==0 & sum(is.na(selCol)==0)) {rownames(df)=unlist(ProteinGeneName[sel, selCol])
  } else {cat("gene lables not updated",sum(is.na(sel)), sum(is.na(selCol)), "\n")}
  #match sampleid order
  new_order=match(colnames(df), annotation$sampleid)
  annotation=annotation[new_order, ]
  data.in <- df
  gene_annot_info <- NULL
  sample_annot=NULL #column annotation
  if (!is.null(heatmap_annot)) {
    sel_col=match(heatmap_annot, names(annotation))
    df_annot=annotation[, sel_col, drop=FALSE]
    sample_annot=HeatmapAnnotation(df = df_annot)
  }
  cluster_rows = FALSE;cluster_cols=FALSE
  if (dendrogram == "both" | dendrogram == "row")
    cluster_rows = TRUE
  if (dendrogram == "both" | dendrogram == "column")
    cluster_cols = TRUE
  
  cexRow = as.numeric(as.character(hyfontsizep))
  cexCol = as.numeric(as.character(hxfontsizep))
  
  labCol = TRUE
  labRow = TRUE
  # cat("pheatmap ", dim(data.in), date(), "\n") #debug
  if (cexRow  == 0 | nrow(data.in) > heatmap_N_genes) {
    labRow = FALSE
    cexRow = 5
  }
  if (cexCol == 0) {
    labCol = FALSE
    cexCol  = 5
  }
  
  cutree_rows = cutreerows
  cutree_cols = cutreecols
  
  #clean up SD=0 rows and columns
  if (scale=="row") {
    row_SD=apply(data.in, 1, function(x) sd(x,na.rm=T))
    data.in=data.in[row_SD!=0, ]
  }
  if (scale=="column") {
    col_SD=apply(data.in, 2, function(x) sd(x,na.rm=T))
    data.in=data.in[, col_SD!=0]
  }	
  
  #now reproduce in Heatmap
  if (scale=="none") {
    data_range=quantile(unlist(data.in), probs=c(0.01, 0.5, 0.99), na.rm=T)
    col_fun=colorRamp2(data_range, c(lowColor,midColor, highColor) )
    legend_text=exp_unit
  } else {
    if (scale=="row") {
      data.in=t(scale(t(data.in)) ) } else {data.in=scale(data.in) }
    data_range=quantile(unlist(abs(data.in)), probs=c(0.01, 0.5, 0.99), na.rm=T)
    max_s=data_range[3]
    col_fun=colorRamp2(c(0-max_s, 0, max_s),  c(lowColor,midColor, highColor) )
    legend_text=str_c("Z-Score of\n", exp_unit)	
  }
  if (cluster_cols==F) {cutree_cols=0}
 # if (heatmap_highlight=="No") {row_label_side="right"} else (row_label_side="left")
  row_label_side="right"
  
  #browser() #debug
  
  p<-Heatmap(data.in, col=col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_cols, 
             clustering_distance_rows=distanceMethod, clustering_distance_columns=distanceMethod,
             clustering_method_rows=agglomerationMethod, clustering_method_columns=agglomerationMethod,
             row_km=cutree_rows, column_km=cutree_cols, row_km_repeats = 100, column_km_repeats = 100,
             show_row_names = labRow, show_column_names = labCol, row_names_side=row_label_side,
             show_row_dend=as.logical(heatmap_row_dend), show_column_dend = as.logical(heatmap_col_dend),
             top_annotation = sample_annot,	row_names_gp = gpar(fontsize = cexRow),
             row_title=heatmap_row_title, column_title=heatmap_column_title,
             row_title_gp = gpar(fontsize = heatmap_row_title_font_size), column_title_gp = gpar(fontsize = heatmap_column_title_font_size),
             column_names_gp = gpar(fontsize = cexCol), heatmap_legend_param = list(title = legend_text, color_bar = "continuous") )
  
  #highlight genes
  if (!is.null(gene_annot_info)) {
    df=gene_annot_info[, 3:ncol(gene_annot_info), drop=F]
    sel_col_path=match(c("Color", "Pathways"), names(df))
    if (sum(is.na(sel_col_path))==0) {
      df_color<-df%>%filter(!duplicated(Color))
      pathway_color=df_color$Color
      names(pathway_color)=df_color$Pathways
      rowAnnot=rowAnnotation(Pathways=gene_annot_info$Pathways, col=list(Pathways=pathway_color) ) 
    } else {rowAnnot=rowAnnotation(df=df)}
    p<-p+rowAnnot
  } 

  return(p)

}
                       
                      
