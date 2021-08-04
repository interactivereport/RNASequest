#!/usr/bin/env Rscript
suppressMessages(require(tidyverse))
suppressMessages(require(reshape2))
suppressMessages(require(DESeq2))
suppressMessages(require(edgeR))
suppressMessages(require(limma))
suppressMessages(require(stringr))
#suppressMessages(require(future.apply))
suppressMessages(require(BiocParallel))

read_file <- function(file, rowname) {
  if (rowname) {
    if (str_detect(file, ".csv$")) {
      table <- read.table(file,header=T,sep=",",as.is=T, check.names=F, row.names =1)
    } else {
      table <- read.table(file,header=T,sep="\t",as.is=T, check.names=F, row.names =1)
    }
  } else {
    if (str_detect(file, ".csv$")) {
      table <- read.table(file,header=T,sep=",",as.is=T, check.names=F)
    } else {
      table <- read.table(file,header=T,sep="\t",as.is=T, check.names=F)
    }
  }
  table[is.na(table)] <- ""
  return(table)
}

Batch_DEG = function(Counts_table, S_meta, comp_info, create_beta_coef_matrix=F,core=2) {
  register(MulticoreParam(core))
  comp_info.list <- setNames(split(cbind(CompareName=rownames(comp_info),comp_info),
                                   seq(nrow(comp_info))),
                             rownames(comp_info))
  #DEG_result_list = list()
  #for (i in names(comp_info.list)) {
  #  DEG_result_list[[i]] = DEG_analysis(comp_info.list[[i]],Counts_table,S_meta,create_beta_coef_matrix)
  #}
  DEG_result_list <- sapply(comp_info.list, DEG_analysis, Counts_table, S_meta, create_beta_coef_matrix)
  return(DEG_result_list)
}


DEG_analysis = function(comp_info,Counts_table,S_meta, create_beta_coef_matrix) {
  comp_name = comp_info$CompareName
  model = comp_info$Model
  group_var = comp_info$Group_name
  subsetting_cov = comp_info$Covariate_subsetting_group
  trt_group = comp_info$Group_test
  ctrl_group = comp_info$Group_ctrl
  analysis_method = comp_info$Analysis_method
  shrink_logFC = comp_info$Shrink_logFC
  LFC_cutoff = comp_info$LFC_cutoff

  S_meta_sub = S_meta
  Counts_table_sub = Counts_table
  
  Subset_group = comp_info$Subsetting_group
  if (!is.na(Subset_group) && !Subset_group == "") {
    Subset = subset_data(Subset_group, S_meta_sub, Counts_table_sub)
    S_meta_sub = Subset$S_meta
    Counts_table_sub =Subset$Counts_table
  }
    
  beta_coef_matrix = data.frame()
  # get covariate factors and adjust basal levels for all the covariates in the model
  if (analysis_method == "DESeq2") {
    for (n in 1: ncol(S_meta_sub)) {
      if(is.numeric(S_meta_sub[,n])) next # O'Young: possible numeric (co-)variates
      S_meta_sub[,n]= factor(S_meta_sub[,n])
    }

    S_meta_sub[,group_var] = relevel(S_meta_sub[,group_var], ref = ctrl_group)
    
    if (!comp_info$Covariate_levels =="") {
      Covariate_levels = strsplit(comp_info$Covariate_levels, "; |;")[[1]]
      Covariate_level_vec = character(length(Covariate_levels))
      Covariate_level_vec_name = character(length(Covariate_levels))
      for (i in 1:length(Covariate_levels)) {
        Covariate_level_vec[i] = strsplit(Covariate_levels[i], ":")[[1]][2]
        Covariate_level_vec_name[i] = strsplit(Covariate_levels[i], ":")[[1]][1]
      }
      names(Covariate_level_vec) = Covariate_level_vec_name
      for (n in 1: length(Covariate_level_vec)) {
        S_meta_sub[,names(Covariate_level_vec)[n]] = relevel(S_meta_sub[,names(Covariate_level_vec)[n]], ref = Covariate_level_vec[n])
      }
    }
    var_list= unique(strsplit(model, "~|\\+| |\\*|\\:")[[1]])
    var_list = var_list[!var_list ==""]
    for (n in 1: length(var_list)) {
      assign(var_list[n], S_meta_sub[, var_list[n]])
    }
    # Design
    dds <- DESeqDataSetFromMatrix(countData = round(Counts_table_sub),
                                  colData = S_meta_sub,#[,unique(trimws(unlist(strsplit(model,"\\+|\\*|\\:"))))]
                                  design= as.formula(paste0('~', model)))
    dds <- DESeq(dds,parallel=T)
    if (LFC_cutoff == 0) {
      res <- results(dds, name = str_c(group_var, "_", trt_group, "_vs_", ctrl_group))
    } else if (LFC_cutoff > 0) {
      res <- results(dds, name = str_c(group_var, "_", trt_group, "_vs_", ctrl_group), lfcThreshold=LFC_cutoff, altHypothesis="greaterAbs")
    }
    # to shrink log fold changes association with condition:
    if (toupper(shrink_logFC) == "YES") {
      res <- lfcShrink(dds, coef=str_c(group_var, "_", trt_group, "_vs_", ctrl_group), type="apeglm",parallel=T)
    }
    if (create_beta_coef_matrix) {
      beta_coef_matrix = beta.coef(assay(vst(dds)), model.matrix(design(dds)), coef(dds), digits = 5)
      # write.csv(beta_coef_matrix, str_c(comp_name, "_beta_coef.csv"))
    }
    
    comp_result = as.data.frame(res[, c(2,4,5)])
    colnames(comp_result) =  paste0(comp_name,"_DESeq.",c("logFC","P.value","Adj.P.value"))
    #rownames(comp_result) = str_split(rownames(comp_result), '\\.', simplify = T)[,1]
    # write.csv(comp_result, str_c(comp_name, "_DEG.csv"))

  } else if (analysis_method == "limma") {
    var_list= unique(strsplit(model, "~|\\+| |\\*|\\:")[[1]])
    var_list = var_list[!var_list ==""]
    for (n in 1: length(var_list)) {
      assign(var_list[n], as.factor(S_meta_sub[, var_list[n]]))
      if (var_list[n] == group_var) {
        assign(var_list[n], relevel(eval(as.name(group_var)), ref = ctrl_group))
      }
    }
    if (!comp_info$Covariate_levels =="") {
      Covariate_levels = strsplit(comp_info$Covariate_levels, "; |;")[[1]]
      Covariate_level_vec = character(length(Covariate_levels))
      Covariate_level_vec_name = character(length(Covariate_levels))
      for (i in 1:length(Covariate_levels)) {
        Covariate_level_vec[i] = strsplit(Covariate_levels[i], ":")[[1]][2]
        Covariate_level_vec_name[i] = strsplit(Covariate_levels[i], ":")[[1]][1]
      }
      names(Covariate_level_vec) = Covariate_level_vec_name
      for (n in 1: length(Covariate_level_vec)) {
        assign(names(Covariate_level_vec)[n], relevel(eval(as.name(names(Covariate_level_vec)[n])), ref = Covariate_level_vec[n]))
      }
    }
    design = model.matrix(as.formula(str_c('~', model)))
    rownames(design) = rownames(S_meta_sub)
    colnames(design) = str_replace(colnames(design), group_var, "")
    column = which(colnames(design) == trt_group)
    
    d2=DGEList(counts = Counts_table_sub)
    d2=calcNormFactors(d2)
    y = voom(d2, design)
    fit = lmFit(y, design)
    
    if (LFC_cutoff == 0) {
      fit2=eBayes(fit)
      tt=topTable(fit2, coef=column, adjust="BH", n=nrow(fit2))
    } else if (LFC_cutoff > 0) {
      fit2<- treat(fit, lfc=LFC_cutoff, trend=TRUE, robust=TRUE)
      tt <- topTreat(fit2, coef=column, adjust.method="BH", n=nrow(fit2))
    }
    if (create_beta_coef_matrix) {
      beta_coef_matrix = beta.coef(y$E, fit2$design, fit2$coefficients, digits = 5)
      # write.csv(beta_coef_matrix, str_c(comp_name, "_beta_coef.csv"))
    }
    if (toupper(shrink_logFC) == "YES") {
      # output predictive log fold changes for first 5 genes
      tt$logFC = predFCm(fit2,coef=column, all.de=T, prop.true.null.method="lfdr")
    }
    ###########################
    #                               logFC    AveExpr         t      P.Value    adj.P.Val        B
    # ENSMUSG00000026822.14 -1.368611e-08 -2.3752664 -8.347657 6.698204e-13 9.400260e-09 17.67702
    # ENSMUSG00000038357.10  1.238756e-08 -1.9859620 -7.518133 3.541754e-11 2.485249e-07 14.32863
    # ENSMUSG00000052212.6  -2.041382e-09 -2.2488520 -7.146580 2.039896e-10 9.542632e-07 12.68047
    
    
    # topTreat table doesn't have B column
    ###########################    
    comp_result = tt[, c(1,4,5)]
    colnames(comp_result) =  paste0(comp_name,"_limma.",c("logFC","P.value","Adj.P.value"))
    #rownames(comp_result) = str_split(rownames(comp_result), '\\.', simplify = T)[,1]
    # write.csv(comp_result, str_c(comp_name, "_DEG.csv"))
  }
  if (nrow(beta_coef_matrix) > 0) {
    result_list <- list("DEG" = comp_result, "beta_coef_matrix" = beta_coef_matrix)
  } else {
    result_list <- list(DEG = comp_result,name=comp_name)
  }
  return(list(result_list))
}

subset_data <- function(Subset_group, Sample_meta, Counts_table) {
  Subset_group_levels = strsplit(Subset_group, "; |;")[[1]]
  Subset_group_level_vec = character(length(Subset_group_levels))
  Subset_group_level_vec_name = character(length(Subset_group_levels))
  for (i in 1:length(Subset_group_levels)) {
    Subset_group_level_vec[i] = strsplit(Subset_group_levels[i], ":")[[1]][2]
    Subset_group_level_vec_name[i] = strsplit(Subset_group_levels[i], ":")[[1]][1]
  }
  names(Subset_group_level_vec) = Subset_group_level_vec_name
  
  for (m in length(Subset_group_level_vec)) {
    Sample_meta = Sample_meta %>% dplyr::filter(get(names(Subset_group_level_vec)[m]) == Subset_group_level_vec[m])
  }
  Counts_table = Counts_table[,rownames(Sample_meta)]
  subset_list = list("S_meta" = Sample_meta, "Counts_table" = Counts_table)
  return(subset_list)
}


beta.coef <- function(response_table, design_table, coef_table, digits = 5) { # start function code
    # response_table, table of response variable, y$E from limma voom() object
    # model  = model from lmfit(), eBayes() or treat()
    # digits = the decimal places to display (default = 5)
    # Code dervied from https://www.dataanalytics.org.uk/wp-content/uploads/2019/08/Beta-coeff-calc.r
    
    bc = data.frame()
    for (i in 1:nrow(response_table)) {
        data = cbind(t(t(response_table[i,])), design_table)
        sdev <- apply(data, 2, sd)
        nvar <- length(sdev)
        for (j in 2:nvar) {
            bc[i,j-1] = coef_table[i,j-1]*sdev[j]/sdev[1]
        }
    }
    rownames(bc) = rownames(coef_table)
    colnames(bc) = colnames(coef_table)
    return(bc)
}


#result = Batch_DEG(rawC, S_meta, comp_info, create_beta_coef_matrix)
