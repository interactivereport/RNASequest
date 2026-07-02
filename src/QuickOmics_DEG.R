#!/usr/bin/env Rscript
suppressMessages(require(tidyverse))
suppressMessages(require(reshape2))
suppressMessages(require(DESeq2))
suppressMessages(require(edgeR))
suppressMessages(require(limma))
suppressMessages(require(stringr))
suppressMessages(require(BiocParallel))

read_file <- function(strF,flag) {
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
  if (grepl("csv$",strF, ".csv$")) {
    comp_info <- as.data.frame(read.table(strF,header=T,sep=",",as.is=T, check.names=F))
  } else {
    comp_info <- as.data.frame(read.table(strF,header=T,sep="\t",as.is=T, check.names=F))
  }
  if(sum(!comTitle%in%colnames(comp_info))>0){
    message("Missing following columns in comparison definition file: ",strF)
    message("\t",paste(comTitle[!comTitle%in%colnames(comp_info)],collapse=", "))
    q()
  }
  comp_info[is.na(comp_info)] <- ""
  comp_info <- data.frame(row.names=unlist(comp_info$CompareName),
                          comp_info[,comTitle[-1]])
  return(comp_info)
}

group_DEG <- function(comp_info){
	# group the compare information by model, subset, Covariate_levels and Analysis_method
	selCol <- c("Model","Subsetting_group","Covariate_levels","Analysis_method")
	grpCom <- comp_info[,selCol]
	grpCom <- grpCom[!duplicated(grpCom),]
	rownames(grpCom) <- paste0("comparisonModel",1:nrow(grpCom))
	grpCom <- apply(grpCom,1,function(x){
		sel <- rep(T,nrow(comp_info))
		for(one in names(x)){
			if(is.na(x[one])) sel <- sel & is.na(unlist(comp_info[,one]))
			else sel <- sel & unlist(comp_info[,one]==x[one])
		}
		return(comp_info[sel,,drop=F])
	})
	names(grpCom) <- rep("",length(grpCom))
	return(grpCom)
}

Batch_DEG = function(Counts_table, S_meta, comp_info, create_beta_coef_matrix=F,core=2) {
  register(MulticoreParam(core))
  #comp_info.list <- setNames(split(cbind(CompareName=rownames(comp_info),comp_info),
  #                                 seq(nrow(comp_info))),
  #                           rownames(comp_info))
  comp_info.list <- group_DEG(cbind(CompareName=rownames(comp_info),comp_info))
  DEG_result_list <- sapply(comp_info.list, DEG_analysis, Counts_table, S_meta, create_beta_coef_matrix)
  DEG_result_list <- unlist(DEG_result_list,recursive=F)
  return(DEG_result_list)
}

DEG_analysis = function(comp_info,Counts_table,S_meta, create_beta_coef_matrix) {
  comp_name = comp_info$CompareName
  Subset_group = unique(comp_info$Subsetting_group)
  analysis_method = unique(comp_info$Analysis_method)
  
  if(length(Subset_group)>1 || length(analysis_method)>1)
  	stop("Error in DEG_analysis comparison group, more than one subset!")
  
  S_meta_sub = S_meta
  Counts_table_sub = Counts_table
  message(paste0("Comparison for ",paste(comp_name,collapse="; ")))
  if (!is.na(Subset_group) && !Subset_group == "") {
    Subset = subset_data(Subset_group, S_meta_sub, Counts_table_sub)
    S_meta_sub = Subset$S_meta
    Counts_table_sub =Subset$Counts_table
  }
  
  # get covariate factors and adjust basal levels for all the covariates in the model
  if (analysis_method == "DESeq2") {
    print(system.time({result_list = DESeq2_DEG(S_meta_sub, Counts_table_sub, comp_info, create_beta_coef_matrix)}))
  } else if (analysis_method == "limma") {
  	#result_list <- sapply(setNames(split(comp_info,seq(nrow(comp_info))),rownames(comp_info)),
  	#					  function(one)return(limma_DEG(S_meta_sub, Counts_table_sub, one, create_beta_coef_matrix)))
  	result_list <- bplapply(setNames(split(comp_info,seq(nrow(comp_info))),rownames(comp_info)),
  							function(one)return(limma_DEG(S_meta_sub, Counts_table_sub, one, create_beta_coef_matrix)))
    #result_list = limma_DEG(S_meta_sub, Counts_table_sub, comp_info, create_beta_coef_matrix)
  	result_list <- unlist(result_list,recursive=F)
  }
  return(list(result_list))
}

subset_data <- function(Subset_group, Sample_meta, Counts_table=NULL) {
    matches <- gregexpr("[^;|]+|[;|]", Subset_group, perl=TRUE)
    subgrp <- trimws(regmatches(Subset_group, matches)[[1]])
    selS <- subset_one(Sample_meta,subgrp[1])
    i<-2
    while(i<length(subgrp)){
        if(subgrp[i]==";"){
            selS <- selS & subset_one(Sample_meta,subgrp[i+1])
        }else if(subgrp[i]=="|"){
            selS <- selS | subset_one(Sample_meta,subgrp[i+1])
        }else{
            stop("unknown seperator: ",subgrp[i])
        }
        i <- i+2
    }
    Sample_meta <- Sample_meta[selS,]
    #Counts_table <- NULL
    if(!is.null(Counts_table))
        Counts_table = Counts_table[,rownames(Sample_meta)]
    return(list(S_meta=Sample_meta,Counts_table=Counts_table))
}
subset_one <- function(sInfo,kv){
    kv = trimws(strsplit(kv,":")[[1]])
    return(sInfo[[kv[1]]]==kv[2])
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

DESeq2_DEG <- function(S_meta, Counts_table, comp_info, create_beta_coef_matrix) {
  comp_info <- comp_info[order(comp_info$Group_ctrl),]
  comp_name = comp_info$CompareName
  model = unique(comp_info$Model)
  
  group_var = unique(comp_info$Group_name)
  trt_group = comp_info$Group_test
  ctrl_group = comp_info$Group_ctrl
  if(nchar(trt_group)==0 || nchar(ctrl_group)==0){
      message("\tNumeric valurable contrast!")
      S_meta[,group_var]= as.numeric(S_meta[,group_var])
  }else{
      S_meta[,group_var]= as.factor(S_meta[,group_var])
  }
  
  Subset_group = unique(comp_info$Subsetting_group)
  Covariate_levels = unique(comp_info$Covariate_levels)
  analysis_method = unique(comp_info$Analysis_method)
  shrink_logFC = comp_info$Shrink_logFC
  LFC_cutoff = comp_info$LFC_cutoff
  
  if(length(model)>1 || length(group_var)>1 || length(Subset_group)>1 || length(Covariate_levels)>1 || length(analysis_method)>1)
  	stop("Error in DESeq2_DEG, more than one comparison subset!")
  
  for (n in colnames(S_meta)) {
    if(n==group_var || is.numeric(S_meta[,n])) next # O'Young: possible numeric (co-)variates
    S_meta[,n]= factor(S_meta[,n])
  }
  #S_meta[,group_var] = relevel(S_meta[,group_var], ref = ctrl_group)
  
  if (str_detect(model, "\\*|\\:") && !Covariate_levels =="") {
    Cov_levels = trimws(strsplit(Covariate_levels, ";")[[1]])
    Cov_levels_vec = setNames(trimws(sapply(strsplit(Cov_levels,":"),tail,1)),trimws(sapply(strsplit(Cov_levels,":"),head,1)))
    for (n in 1: length(Cov_levels_vec)) {
      S_meta[,names(Cov_levels_vec)[n]] = relevel(S_meta[,names(Cov_levels_vec)[n]], ref = Cov_levels_vec[n])
    }
  }
  #var_list= unique(as.vector(trimws(str_split(model, "~|\\+|\\*|\\:", simplify =T))))
  #var_list = var_list[!var_list ==""]
  #for (n in 1: length(var_list)) {
  #  assign(var_list[n], S_meta[, var_list[n]])
  #}
  # Design
  dds <- DESeqDataSetFromMatrix(countData = round(Counts_table),
                                colData = S_meta,#[,unique(trimws(unlist(strsplit(model,"\\+|\\*|\\:"))))]
                                design= as.formula(paste0('~', model)))
  dds <- DESeq(dds,parallel=T)
  result_list <- list()
  
  for(i in 1:nrow(comp_info)){
  	message("\t--- ",comp_name[i])#
    if(is.numeric(S_meta[,group_var]))
        res <- results(dds,name=group_var,
                       lfcThreshold=LFC_cutoff[i], altHypothesis="greaterAbs",
                       parallel=T)
    else
        res <- results(dds,contrast=c(group_var,trt_group[i],ctrl_group[i]),
                       lfcThreshold=LFC_cutoff[i], altHypothesis="greaterAbs",
                       parallel=T)
    
  	if (toupper(shrink_logFC[i]) == "YES") {
  	    if(is.numeric(S_meta[,group_var])) strContrast <- group_var
  	    else strContrast <- str_c(group_var, "_", trt_group[i], "_vs_", ctrl_group[i])
  		if(!strContrast%in%resultsNames(dds)){
  			#https://www.biostars.org/p/448959/
  			dds[[group_var]] <- relevel(dds[[group_var]], ctrl_group[i])
  			dds <- nbinomWaldTest(dds)
  			message("\t\trelevel")
  		}
  		shrinkLFC <- suppressMessages(lfcShrink(dds,coef=strContrast, type="apeglm",
  		                                  format="DataFrame",svalue=T,
  		                                  parallel=T,lfcThreshold=LFC_cutoff[i]))
  		lfc <- "log2FoldChange"
  		res[[lfc]] <- shrinkLFC[rownames(res),lfc]
  		res <- cbind(res,data.frame(svalue=shrinkLFC[rownames(res),"svalue"],
  		                            sadj=p.adjust(shrinkLFC[rownames(res),"svalue"],"BH")))
  	}
  	if (create_beta_coef_matrix) {
  		beta_coef_matrix = data.frame()
  		beta_coef_matrix = beta.coef(assay(vst(dds)), model.matrix(design(dds)), coef(dds), digits = 5)
  		# write.csv(beta_coef_matrix, str_c(comp_name, "_beta_coef.csv"))
  	}
  	sel <- c(grep("log",colnames(res),ignore.case=T,value=T),
  	         grep("pvalue",colnames(res),ignore.case=T,value=T),
  	         grep("padj",colnames(res),ignore.case=T,value=T),
  	         grep("svalue",colnames(res),ignore.case=T,value=T),
  	         grep("sadj",colnames(res),ignore.case=T,value=T))
  	comp_result = as.data.frame(res[, sel])#log2FoldChange, pvalue, padj
  	colnames(comp_result) =  paste0(comp_name[i],"_DESeq2.",colnames(comp_result))
  	if (exists("beta_coef_matrix") && nrow(beta_coef_matrix) > 0) {
  		result_list[[length(result_list)+1]] <- list("DEG" = comp_result, "beta_coef_matrix" = beta_coef_matrix)
  	} else {
  		result_list[[length(result_list)+1]] <- list(DEG = comp_result,name=comp_name[i])
  	}
  }
  names(result_list) <- rownames(comp_info)
  return(result_list)
  
  # if (LFC_cutoff == 0) {
  #   res <- results(dds, name = str_c(group_var, "_", trt_group, "_vs_", ctrl_group))
  # } else if (LFC_cutoff > 0) {
  #   res <- results(dds, name = str_c(group_var, "_", trt_group, "_vs_", ctrl_group), lfcThreshold=LFC_cutoff, altHypothesis="greaterAbs")
  # }
  # 
  # # to shrink log fold changes association with condition:
  # if (toupper(shrink_logFC) == "YES") {
  #   res <- lfcShrink(dds, coef=str_c(group_var, "_", trt_group, "_vs_", ctrl_group), type="apeglm", parallel=T)
  # }
  # if (create_beta_coef_matrix) {
  #   beta_coef_matrix = data.frame()
  #   beta_coef_matrix = beta.coef(assay(vst(dds)), model.matrix(design(dds)), coef(dds), digits = 5)
  #   # write.csv(beta_coef_matrix, str_c(comp_name, "_beta_coef.csv"))
  # }
  # 
  # comp_result = as.data.frame(res[, c(2,4,5)])
  # colnames(comp_result) =  paste0(comp_name,"_DESeq.",c("logFC","P.value","Adj.P.value"))
  # #rownames(comp_result) = str_split(rownames(comp_result), '\\.', simplify = T)[,1]
  # # write.csv(comp_result, str_c(comp_name, "_DEG.csv"))
  # 
  # if (exists("beta_coef_matrix") && nrow(beta_coef_matrix) > 0) {
  #   result_list <- list("DEG" = comp_result, "beta_coef_matrix" = beta_coef_matrix)
  # } else {
  #   result_list <- list(DEG = comp_result,name=comp_name)
  # }
  # return(list(result_list))
} 

limma_DEG <- function(S_meta, Counts_table, comp_info, create_beta_coef_matrix) {
  comp_name = comp_info$CompareName
  message("\t--- ",comp_name)
  model = comp_info$Model
  group_var = comp_info$Group_name
  Subset_group = comp_info$Subsetting_group
  Covariate_levels = comp_info$Covariate_levels
  trt_group = comp_info$Group_test
  ctrl_group = comp_info$Group_ctrl
  analysis_method = comp_info$Analysis_method
  #shrink_logFC = comp_info$Shrink_logFC
  LFC_cutoff = comp_info$LFC_cutoff
  
  var_list= unique(trimws(str_split(model, "~|\\+|\\*|\\:", simplify =T)[1,]))
  var_list = var_list[!var_list ==""]
  for (n in 1: length(var_list)) {
    if(is.numeric(S_meta[,var_list[n]])) {
      assign(var_list[n], S_meta[, var_list[n]])
    } else {
      assign(var_list[n], as.factor(S_meta[, var_list[n]]))
    }
  }
  
  if (str_detect(model, "\\*|\\:") && !Covariate_levels =="") {
    Cov_levels = trimws(strsplit(Covariate_levels, ";")[[1]])
    Cov_levels_vec = setNames(trimws(sapply(strsplit(Cov_levels,":"),tail,1)),trimws(sapply(strsplit(Cov_levels,":"),head,1)))
    for (n in 1: length(Cov_levels_vec)) {
      assign(names(Cov_levels_vec)[n], relevel(eval(as.name(names(Cov_levels_vec)[n])), ref = Cov_levels_vec[n]))
    }
  }

  assign(group_var, relevel(eval(as.name(group_var)), ref = ctrl_group))
  
  design = model.matrix(as.formula(str_c('~', model)))
  rownames(design) = rownames(S_meta)
  colnames(design) = str_replace(colnames(design), group_var, "")
  column = which(colnames(design) == trt_group)
  
  d2=DGEList(counts = Counts_table)
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
    beta_coef_matrix = data.frame()
    beta_coef_matrix = beta.coef(y$E, fit2$design, fit2$coefficients, digits = 5)
    # write.csv(beta_coef_matrix, str_c(comp_name, "_beta_coef.csv"))
  }
  #if (toupper(shrink_logFC) == "YES") {
  #  # output predictive log fold changes for first 5 genes
  #  tt$logFC = predFCm(fit2,coef=column, all.de=T, prop.true.null.method="lfdr")
  #}
  
  ###########################  
  sel <- c(grep("log",colnames(tt),ignore.case=T,value=T),
           grep("^p.",colnames(tt),ignore.case=T,value=T),
           grep("^adj.",colnames(tt),ignore.case=T,value=T))
  comp_result = tt[, sel]#"logFC","P.value","Adj.P.value"
  colnames(comp_result) =  paste0(comp_name,"_limma.",c("logFC","P.value","Adj.P.value"))
  #rownames(comp_result) = str_split(rownames(comp_result), '\\.', simplify = T)[,1]
  # write.csv(comp_result, str_c(comp_name, "_DEG.csv"))
  
  if (exists("beta_coef_matrix") && nrow(beta_coef_matrix) > 0) {
    result_list <- list("DEG" = comp_result, "beta_coef_matrix" = beta_coef_matrix)
  } else {
    result_list <- list(DEG = comp_result,name=comp_name)
  }
  return(list(result_list))
}

checkComparisonInfo <- function(comp_info, meta, comp_info_file) {
  requiredCols = c("Subsetting_group","Model","Covariate_levels","Group_name",
				   "Group_test","Group_ctrl","Analysis_method","Shrink_logFC","LFC_cutoff")
  notmatching = setdiff(requiredCols,names(comp_info))
  if (length(notmatching) > 0) 
  { 
    stop(paste0(comp_info_file," should have column header(s): ",
         paste(notmatching, collapse=',')
		 ))
  }
  if(nrow(comp_info)==0){
    stop(paste0("Empty comparison definition file (", comp_info_file, ") is NOT allowed!"))
  }
  
  # set the default value if those are empty
  setDefault <- F
  for(i in rownames(comp_info)){
    if(is.null(comp_info[i,"Group_name"]) || nchar(comp_info[i,"Group_name"])==0){
      stop(paste0("'Group_name' cannot be empty for comparison",i))
    }
    if(is.null(comp_info[i,"Model"]) || nchar(comp_info[i,"Model"])==0){
      comp_info[i,"Model"] <- comp_info[i,"Group_name"]
      setDefault <- T
    }
    if(!grepl(comp_info[i,"Group_name"],comp_info[i,"Model"])){
      stop(paste("The DEG 'Model' didn't include 'Group_name' in the comparison file for",i))
    }
    if(is.null(comp_info[i,"Shrink_logFC"]) || nchar(comp_info[i,"Shrink_logFC"])==0){
      comp_info[i,"Shrink_logFC"] <- "Yes"
      setDefault <- T
    }
    if(is.null(comp_info[i,"LFC_cutoff"]) || nchar(comp_info[i,"LFC_cutoff"])==0){
      comp_info[i,"LFC_cutoff"] <- 0
      setDefault <- T
    }
    if(!comp_info[i,"Analysis_method"]%in%c("DESeq2","limma")){
      stop(paste(comp_info[i,"Analysis_method"],"for comparison",i,
                 "is NOT a valide 'Analysis_method' (DESeq2 or limma) in comparison file."))
    }
    if((nchar(comp_info[i,"Group_test"])>0) != (nchar(comp_info[i,"Group_ctrl"])>0)){
        stop(paste("Missing either 'Group_test' or 'Group_ctrl', leave both empty for the numeric 'Group_name'"))
    }
    if(!grepl("^[A-Zza-z0-9._]+$", comp_info[i,"Group_test"]) && nchar(comp_info[i,"Group_test"])>0){
      stop(paste("'Group_test' entry", comp_info[i,"Group_test"],"for comparison",i,
                 "contains characters other than letters, numbers, and delimiters '_' or '.'. Please use only letters, numbers, '_' or '.', as these are safe characters for column names in R."))
    }
    if(!grepl("^[A-Zza-z0-9._]+$", comp_info[i,"Group_ctrl"]) && nchar(comp_info[i,"Group_ctrl"])>0){
      stop(paste("'Group_ctrl' entry", comp_info[i,"Group_ctrl"],"for comparison",i,
                 "contains characters other than letters, numbers, and delimiters '_' or '.'. Please use only letters, numbers, '_' or '.', as these are safe characters for column names in R."))
    }
  }
  comp_info$LFC_cutoff <- as.numeric(comp_info$LFC_cutoff)
  if(sum(comp_info$LFC_cutoff<0)>0) stop("'LFC_cutoff' in comparison file is required to be non-negative!")
  comp_info$Group_test <- as.character(comp_info$Group_test)
  comp_info$Group_ctrl <- as.character(comp_info$Group_ctrl)
  if(setDefault){
    A <- cbind(CompareName=rownames(comp_info),comp_info)
    file.rename(comp_info_file, paste0(comp_info_file,".bk"))
    write.csv(A,file=comp_info_file,row.names=F)
    message("-----> comparison file (",basename(comp_info_file),") is updated with some default values!")
    message("\tThe original comparison file is renamed as ...bk")
  }
  checkComparisonModel(comp_info, meta)
  return(comp_info)
}

checkComparisonModel <- function(comp_info, meta) {
  for (i in 1: nrow(comp_info)) {
    comp_name = rownames(comp_info)[i]
    model = comp_info$Model[i]
    group_var = comp_info$Group_name[i]
    Subset_group = comp_info$Subsetting_group[i]
    Covariate_levels = comp_info$Covariate_levels[i]
    trt_group = comp_info$Group_test[i]
    ctrl_group = comp_info$Group_ctrl[i]
    
    message("\tChecking model for ",comp_name)
    # check the existence of subset variables and subset levels in sample meta table
    if (!is.na(Subset_group) && !Subset_group == "") {
      Subset_group_levels = trimws(strsplit(Subset_group, ";|\\|")[[1]])
      Subset_group_level_vec <- setNames(trimws(sapply(strsplit(Subset_group_levels,":"),tail,1)),trimws(sapply(strsplit(Subset_group_levels,":"),head,1)))
      if(sum(!names(Subset_group_level_vec)%in%colnames(meta))>0) {
        stop(paste("Error in", comp_name, ": Subsetting covariate:",names(Subset_group_level_vec)[!names(Subset_group_level_vec)%in%colnames(meta)], "is NOT defined in the sample meta file"))}
      else {
        Sample_meta = subset_data(Subset_group,meta)$S_meta
        #for(j in names(Subset_group_level_vec)){
        #  if(!Subset_group_level_vec[j]%in%Sample_meta[,j]) {
        #    stop(paste("Error in", comp_name, ": the subsetting covariate value", Subset_group_level_vec[j], "is NOT defined in the", j, "column in the sample meta file"))}
        #  else {
        #    Sample_meta = Sample_meta[Sample_meta[,j]==Subset_group_level_vec[j],]}
        #}
        if(nrow(Sample_meta)==0)
            stop("\n\tNo sample left after subset by: ",Subset_group,"\t(Please note: ';' means 'and', '|' means 'or')\n")
        message("\t\t",nrow(Sample_meta)," samples left after filtering by ",Subset_group)
      }
    } else {Sample_meta =meta}
    
    # check if the comparsion covariate and levels exist in the sample meta table
    if (!group_var %in% names(meta)) {
        stop(paste("Error in", comp_name, ": group variable", group_var, "is NOT defined in the sample meta file."))
    } else if(nchar(trt_group)>0 || nchar(ctrl_group)>0){
        group_var_levels = unique(Sample_meta[,group_var])
        if (!trt_group%in%group_var_levels) {
            stop(paste("Error in", comp_name, "trt_group", trt_group, "is NOT defined in the sample meta file."))
        } else if (!ctrl_group%in%group_var_levels) {
            stop(paste("Error in", comp_name, "ctrl_group", ctrl_group, "is NOT defined in the sample meta file."))
        }
    }
    
    # check if the comparison groups contain multiple samples
    # expect numeric comparison variable
    if(nchar(trt_group)>0 || nchar(ctrl_group)>0){
        n_samples = apply(Sample_meta[, group_var,drop=F], 2, function(x) return(table(x)))
        if (! trt_group %in% rownames(n_samples)) {
            stop(paste("Error in", comp_name, ": trt_group ", trt_group, " is NOT defined in the sample meta file."))
        } else if (n_samples[trt_group,] < 2) {
            stop(paste("Error in", comp_name, ": there is no replicate samples for trt_group ", trt_group, "."))
        }
        
        if (! ctrl_group %in% rownames(n_samples)) {
            stop(paste("Error in", comp_name, ": ctrl_group ", ctrl_group, " is NOT defined in the sample meta file."))
        } else if (n_samples[ctrl_group,] < 2) {
            stop(paste("Error in", comp_name, ": there is no replicate samples for ctrl_group ", ctrl_group, "."))
        }
    }

    # check if all the covariates in the model exist in sample meta table and have multiple levels and are independent
    covariate_list = unique(as.vector(trimws(str_split(model, "~|\\+|\\*|\\:", simplify =T))))
    covariate_list = covariate_list[!covariate_list ==""]
    if (sum(!covariate_list%in%names(meta)) > 0) {
      stop(paste("Error in", comp_name, ": Covariates:", paste(covariate_list[!covariate_list%in%names(meta)], collapse = ','), "is NOT defined in the sample meta file."))
    } else {
      Sample_meta = Sample_meta[, covariate_list,drop=F]
      n_covariate_levels = apply(Sample_meta, 2, function(x) return(length(unique(x))))
      single_level_covariates = names(Sample_meta)[n_covariate_levels == 1]
      if (length(single_level_covariates) > 0) {
        stop(paste("Error in", comp_name, ": Covariates:", paste(single_level_covariates, collapse = ','), "only have one level."))
      } else {
        suppressWarnings(suppressMessages(library(caret)))
        Sample_meta = data.frame(Sample_meta)
        for (n in 1: length(covariate_list)) {
          if(is.numeric(Sample_meta[,covariate_list[n]])) {
            assign(covariate_list[n], Sample_meta[, covariate_list[n]])
          } else {
            assign(covariate_list[n], as.factor(Sample_meta[, covariate_list[n]]))
          }
        }
        
        design = model.matrix(as.formula(str_c('~', model)))
        LinearCombos_result = findLinearCombos(design)
        if (length(LinearCombos_result$linearCombos) > 0) {
          LinearCombos_vec = colnames(design)[LinearCombos_result$linearCombos[[1]]]
          Cov_remove = colnames(design)[LinearCombos_result$remove]
          LinearCombos_list = character()
          Cov_remove_list = character()
          for (var in covariate_list) {
            if (sum(grepl(var, LinearCombos_vec)) > 0 )
              LinearCombos_list = c(LinearCombos_list, var)
            if (sum(grepl(var, Cov_remove)) > 0 )
              Cov_remove_list = c(Cov_remove_list, var)
          }
          #stop(paste("Error in", comp_name, "model:", paste(LinearCombos_list, collapse = ' and '), "are linear_combinated variables. Please remove", paste(Cov_remove_list, collapse = ',')))          
        }
      }
    }
    
    
    # check the existence of varibales in the interaction term and the specified levels in sample meta table
    if (str_detect(model, "\\*|\\:")) {
      # Covariate_levels is required if group_var is involved in interactions with model with interaction terms.
      model_additive = trimws(str_split(model, "\\+", simplify =T)) # separate the additive covariates
      group_var_substr = model_additive[str_detect(model_additive, group_var)]  # get additive model element contains the group_var
      if (sum(str_detect(group_var_substr, "\\*|\\:"))>0) {
        interaction_terms = group_var_substr[str_detect(group_var_substr, "\\*|\\:")]
        interaction_covariates = setdiff(trimws(str_split(interaction_terms, "\\*|\\:", simplify =T)), group_var) # get list of covariates other than group_var in the interaction terms.
        
        var_levels = trimws(strsplit(Covariate_levels, ";")[[1]])
        var_level_vec = setNames(trimws(sapply(strsplit(var_levels,":"),tail,1)),trimws(sapply(strsplit(var_levels,":"),head,1)))
        if (sum(!interaction_covariates%in%names(var_level_vec))>0)  {
          stop(paste("Error in", comp_name, ": covariate_levels are required to set for", paste(interaction_covariates, collapse = ','), ", which has interaction in the model."))
        } else if(sum(!names(var_level_vec)%in%colnames(meta))>0) {
          stop(paste("Error in", comp_name, ": Covariate_levels:", names(var_level_vec)[!names(var_level_vec)%in%colnames(meta)], "is NOT defined in the sample meta file"))
        } else {
          for(j in names(var_level_vec)){
            if(!var_level_vec[j]%in%Sample_meta[,j]) {
              stop(paste("Error in", comp_name, ": Covariate", j, " level", var_level_vec[j], "does not exist in the", names(var_level_vec)[j], "column in the sample meta file"))}
          }      
        }
      }
    }
  } 
}

