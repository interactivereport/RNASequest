#! /usr/bin/env Rscript
options(stringsAsFactors=F)
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(stringr))

rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args)<2){
  stop("The config file of a finished EArun project is required!")
}

## -----
message("loading resource ...")
config <- yaml::read_yaml(args[2])
sys_config <- yaml::read_yaml(file.path(args[1],"sys.yml"))
source(file.path(args[1],"QuickOmics_DEG.R"))
# config <- yaml::read_yaml("/camhpc/ngs/projects/TST11773/dnanexus/test_zhengyu.ouyang/QuickOmics_20210712180059/config.yml")
# config <- yaml::read_yaml("/camhpc/ngs/projects/TST11589/dnanexus/20210426220540_zhengyu.ouyang/QuickOmics_20210708004132/config.yml")
# sys_config <- yaml::read_yaml(file.path("/home/wli7/projects/Quickomics_converter/src","sys.yml"))
# source(file.path("/home/wli7/projects/Quickomics_converter/src","QuickOmics_DEG.R"))

R_file = normalizePath(file.path(sys_config$QuickOmics_path, paste0(config$prj_name, ".RData")))
if(!file.exists(R_file)){
  stop(paste(R_file, "is not found."))
}

project_file = gsub("RData", "csv", R_file)
if(!file.exists(project_file)){
  stop(paste(project_file, "is not found."))
}

## project info -----
p_info = read.csv(project_file)
if (!p_info$Species[1] %in% c("human", "mouse", "rat")) {
  stop("The DA system does not support species other than human, mouse and rat. The program is stopped.")
}

if(is.null(config$DA_file_outpath) || length(config$DA_file_outpath)==0) {
  output_path = file.path(dirname(args[2]), "DA_Import_Files")
} else {
  output_path = config$DA_file_outpath
}
  
# output_path = file.path(dirname("/camhpc/ngs/projects/TST11589/dnanexus/20210426220540_zhengyu.ouyang/QuickOmics_20210708004132/config.yml"), "DA_Import_Files")
dir.create(output_path, showWarnings = FALSE)

project_info  = data.frame(ProjectID=p_info$ProjectID[1], Title=p_info$ShortName[1], Description=p_info$Name[1], Platform='', PlatformDescription='', PlatformTechnology='', ContactName='')
write.csv(project_info, file.path(output_path, "Project_Info.csv"), row.names=F)

## sample info -----
load(R_file)

if (length(which(str_detect(names(MetaData), "Sapio"))) > 0) {
  MetaData = MetaData[, -which(str_detect(names(MetaData), "Sapio"))]
}
col_exclude_list1 = c("Order", "ComparePairs", "Index_ID2", "Annotated_By", "Index_ID", "Status", "Index_Tag", "Index_Tag2", "Status_Time", "Well", "Well_Row", "Plate_Name", "Well_Column")
MetaData = MetaData[, -which(colnames(MetaData) %in% col_exclude_list1)]

# sys_config <- yaml::read_yaml(file.path(args[1],"sys.yml"))
sys_config <- yaml::read_yaml(file.path("/home/wli7/projects/Quickomics_converter/src","sys.yml"))


# col_exclude_list2 = gsub(" ", "_", sys_config$qc2meta)
# MetaData = MetaData[, -which(colnames(MetaData) %in% col_exclude_list2)]

# df_header = read.csv(paste0(args[1],"NGSone2DA_header_mapping.csv"))
df_header = read.csv('/home/wli7/projects/Quickomics_converter/src/NGSone2DA_header_mapping.csv')
df_header = df_header[which(df_header$NGSone_header %in% names(MetaData)),]

meta_names_ori = names(MetaData)
for (i in 1:nrow(df_header)) {
  if (df_header$NGSone_header[i] %in% names(MetaData)) {
    if (!df_header$DA_header[i] %in% names(MetaData)) {
      names(MetaData) = gsub(df_header$NGSone_header[i], df_header$DA_header[i], names(MetaData))
      # MetaData = MetaData %>% dplyr::rename(as.name(df_header$DA_header[i])=as.name(df_header$NGSone_header[i]))
    } else if (!df_header$DA_alt_header[i] == "") {
      names(MetaData) = gsub(df_header$NGSone_header[i], df_header$DA_alt_header[i], names(MetaData))
      # MetaData = MetaData %>% dplyr::rename(df_header$DA_alt_header[i] = df_header$NGSone_header[i])
    }
  }
}

if ("Disease" %in% names(MetaData) && !"DiseaseStage" %in% names(MetaData)) {
  names(MetaData) = gsub("Disease", "DiseaseStage", names(MetaData))
}

DA_header_inclusion_list = sys_config$DA_columns
MetaData = MetaData[, which(names(MetaData) %in% DA_header_inclusion_list)]
MetaData$SampleID = str_replace_all(MetaData$SampleID, "-", "_")
MetaData$SampleID[grepl("^[[:digit:]]+", MetaData$SampleID)] = str_c("S", MetaData$SampleID[grepl("^[[:digit:]]+", MetaData$SampleID)])
write.csv(MetaData, file.path(output_path, "Sample_Info.csv"), row.names=F)

## Expression & Counts ------
if(!is.null(config$prj_path)){
  estCount <- read.table(paste0(config$prj_path,"/combine_rsem_outputs/genes.estcount_table.txt"),
                         header=T,row.names=1,sep="\t",check.names=F,as.is=T)
  estCount <- estCount[!grepl("^ERCC",rownames(estCount)),]
  colnames(estCount) <- sapply(strsplit(sapply(strsplit(colnames(estCount),"\\|"),head,1),
                                        "_"),tail,1)
  if (file.exists(paste0(config$prj_path,"/combine_rsem_outputs/genes.fpkm_table.txt"))) {
    Expression <- read.table(paste0(config$prj_path,"/combine_rsem_outputs/genes.fpkm_table.txt"),
                             header=T,row.names=1,sep="\t",check.names=F,as.is=T)
    
  } else {
    Expression <- read.table(paste0(config$prj_path,"/combine_rsem_outputs/genes.tpm_table.txt"),
                             header=T,row.names=1,sep="\t",check.names=F,as.is=T)
  }
  Expression <- Expression[!grepl("^ERCC",rownames(Expression)),]
  colnames(Expression) <- sapply(strsplit(sapply(strsplit(colnames(Expression),"\\|"),head,1),
                                        "_"),tail,1)
}else{
  if(!is.null(config$exp_counts))
    estCount <- read.table(config$exp_counts,header=T,row.names=1,sep="\t",check.names=F,as.is=T)
  if(!is.null(config$exp_tpm)){
    Expression <-  read.table(config$exp_tpm,header=T,row.names=1, sep="\t",check.names=F,as.is=T)
  }
}
colnames(estCount) = str_replace_all(colnames(estCount), "-", "_")
colnames(Expression) = str_replace_all(colnames(Expression), "-", "_")
colnames(estCount)[grepl("^[[:digit:]]+", colnames(estCount))] = str_c("S", colnames(estCount)[grepl("^[[:digit:]]+", colnames(estCount))])
colnames(Expression)[grepl("^[[:digit:]]+", colnames(Expression))] = str_c("S", colnames(Expression)[grepl("^[[:digit:]]+", colnames(Expression))])


if(!is.null(estCount)) estCount <- estCount[, MetaData$SampleID]
estCount <- estCount[apply(estCount,1,function(x)return(sum(x>=config$min_count)))>=config$min_sample,]

if(!is.null(Expression)) Expression <- Expression[, MetaData$SampleID]

fwrite(data.frame(Gene = rownames(Expression), Expression), file.path(output_path, "Gene_Expression_Data.csv"))
fwrite(data.frame(Gene = rownames(estCount), estCount), file.path(output_path, "Gene_Count.csv"))

## comparison info ------
ComparisonID=rownames(comp_info)

c_info=data.frame(ComparisonID, ProjectName=project_info$ProjectID[1], Case_SampleIDs='', Control_SampleIDs='')
for (i in 1:nrow(c_info)) {
  S_meta_sub = MetaData
  Subset_group = comp_info$Subsetting_group[i]
  if (!is.na(Subset_group) && !Subset_group == "") {
    Subset = subset_data(Subset_group, MetaData, estCount)
    S_meta_sub = Subset$S_meta
  }
  
  if (comp_info$Analysis_method[i] == "DESeq2") {
    c_info$ComparisonID[i] = paste0(c_info$ComparisonID[i], "_DESeq")
  } else if (comp_info$Analysis_method[i] == "limma") {
    c_info$ComparisonID[i] = paste0(c_info$ComparisonID[i], "_limma")
  }    

  compare_group = comp_info$Group_name[i]
  if (!compare_group == "SampleType") {
    c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
    c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
  }

  compare_group_new = compare_group
  if (compare_group == "Age") {
    compare_group_new = "AgeCategory"
  } else if ((compare_group == "Infection" && "Cell_Line" %in% meta_names_ori) || (compare_group == "SampleSource" && !"Organ" %in% meta_names_ori)) {
    compare_group_new = "CellLine"
  } else if (compare_group == "Transfection" && "Genotype" %in% meta_names_ori) {
    compare_group_new = "Genotype"
  } else if (compare_group == "SampleType" && "group" %in% meta_names_ori) {
    compare_group_new = "SubjectGroup"
  }
  c_info[i, paste0("Case_", compare_group_new)] = comp_info$Group_test[i]
  c_info[i, paste0("Control_", compare_group_new)] = comp_info$Group_ctrl[i]
  
  c_info$Case_SampleIDs[i]=str_c(S_meta_sub$SampleID[S_meta_sub[,compare_group]==comp_info$Group_test[i]], collapse=";") 
  c_info$Control_SampleIDs[i]=str_c(S_meta_sub$SampleID[S_meta_sub[,compare_group]==comp_info$Group_ctrl[i]], collapse=";") 
  
  Covariate_list = setdiff(names(S_meta_sub), c("SampleID", "ProjectName", "Organism", compare_group))
  if (length(Covariate_list) > 0) {
    for (j in 1:length(Covariate_list)) {
      case = unique(S_meta_sub[S_meta_sub[,compare_group]==comp_info$Group_test[i], Covariate_list[j]])
      control = unique(S_meta_sub[S_meta_sub[,compare_group]==comp_info$Group_ctrl[i], Covariate_list[j]])

      Covariate = Covariate_list[j]
      if (Covariate_list[j] == "Age") {
        Covariate = "AgeCategory"
      } else if ((Covariate_list[j] == "Infection" && "Cell_Line" %in% meta_names_ori) || (Covariate_list[j] == "SampleSource" && !"Organ" %in% meta_names_ori)) {
        Covariate = "CellLine"
      } else if (Covariate_list[j] == "Transfection" && "Genotype" %in% meta_names_ori) {
        Covariate = "Genotype"
      } else if (Covariate_list[j] == "SampleType" && "group" %in% meta_names_ori) {
        Covariate = "SubjectGroup"
      }
      
      if (length(case)== 1) {
        if (!Covariate_list[j] == "SampleType") {
          c_info[i, paste0("Case_", Covariate_list[j])]=case
        }
        c_info[i, paste0("Case_", Covariate)]=case
      }
      if (length(control)== 1) {
        if (!Covariate_list[j] == "SampleType") {
          c_info[i, paste0("Control_", Covariate_list[j])]=control
        }
        c_info[i, paste0("Control_", Covariate)]=control
      }
    }
  }
} 

fwrite(c_info, file.path(output_path, "Comparisons_Info.csv"))


## comparison data ------
comp_data = data_results[, which(str_detect(names(data_results), "logFC|P.value|Adj.P.value"))]

Comparison_long <- NULL
for (i in seq(1, ncol(comp_data), 3)) {
  subdata=comp_data[, i:(i+2)]
  subdata$Gene=rownames(comp_data)
  subdata$ComparisonID=str_replace(names(subdata)[which(str_detect(names(subdata), "logFC"))], ".logFC", "")
  names(subdata)=str_replace(names(subdata), subdata$ComparisonID[1], "")
  names(subdata)=str_replace(names(subdata), "^\\.", "")
  names(subdata)=str_replace(names(subdata), "Adj.P.value", "adj.P.Val")
  names(subdata)=str_replace(names(subdata), "P.value", "P.Value")
  subdata = subdata[, c('Gene', 'ComparisonID', 'logFC', 'P.Value', 'adj.P.Val')]
  if (is.null(Comparison_long)) {Comparison_long=subdata} else {
    Comparison_long=rbind(Comparison_long, subdata)
  }
}

fwrite(Comparison_long, file.path(output_path, "Comparisons_Data.csv"))


## finished -----
message("=================================================\nResults are saved in ",output_path)
sink(paste0(config$output,"/session.EA2DA"))
sessionInfo()
sink()
