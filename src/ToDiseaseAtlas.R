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

R_file = normalizePath(file.path(config$output, paste0(config$prj_name, ".RData")))
if(!file.exists(R_file)){
  stop(paste(R_file, "is not found."))
}

project_file = gsub("RData", "csv", R_file)
if(!file.exists(project_file)){
  stop(paste(project_file, "is not found."))
}

## project info -----
p_info = read.csv(project_file)
if (!p_info$Species[1] %in% c("human", "homo_sapiens", "mouse", "mus_musculus", "rat", "rattus_norvegicus")) {
  stop("The DA system does not support species other than human, mouse and rat. The program is stopped.")
}

if(is.null(config$DA_file_outpath) || length(config$DA_file_outpath)==0) {
  output_path = file.path(dirname(args[2]), "DA_Import_Files")
} else {
  output_path = config$DA_file_outpath
}
  
dir.create(output_path, showWarnings = FALSE)

project_info  = data.frame(ProjectID=p_info$ProjectID[1], Title=p_info$ShortName[1], Description=p_info$Name[1], ExperimentType="Expression profiling by high throughput sequencing", Platform='', PlatformDescription='', PlatformTechnology='', ContactName='')
write.csv(project_info, file.path(output_path, "Project_Info.csv"), row.names=F)

## sample info -----
load(R_file)
MetaData[] =lapply(MetaData, as.character)
header_mapping = read.csv(file.path(args[1],"NGSone2DA_header_mapping.csv"))
# only keep mapping headers that exist in the MetaData columns and do have mapped DA fields.
header_mapping = header_mapping[intersect(which(header_mapping$NGSone_header %in% names(MetaData)), which(!header_mapping$DA_header == "")),]

# Only keep the default NGSone columns of MetaData that can be mapped to DA
DA_MetaData = MetaData[, which(colnames(MetaData) %in% header_mapping$NGSone_header)]
MetaData_extra = MetaData[, -which(colnames(MetaData) %in% header_mapping$NGSone_header)]
MetaData_extra = MetaData_extra[, -which(colnames(MetaData_extra) %in% c("Order", "ComparePairs"))]

meta_names_ori = names(DA_MetaData)

if (nrow(header_mapping) > 0) {
  for (i in 1:nrow(header_mapping)) {
    if (!header_mapping$NGSone_header[i] == header_mapping$DA_header[i]) {
      if (!header_mapping$DA_header[i] %in% names(DA_MetaData)) {
        names(DA_MetaData) = gsub(header_mapping$NGSone_header[i], header_mapping$DA_header[i], names(DA_MetaData))
      } else {
        DA_MetaData[, header_mapping$DA_header[i]]=str_c(DA_MetaData[, header_mapping$DA_header[i]], "__", DA_MetaData[, header_mapping$NGSone_header[i]])
        DA_MetaData[, header_mapping$DA_header[i]]=str_replace(DA_MetaData[, header_mapping$DA_header[i]], "^__", "")
        DA_MetaData[, header_mapping$DA_header[i]]=str_replace(DA_MetaData[, header_mapping$DA_header[i]], "__$", "")
      }
    }
  }
  if ("DiseaseState" %in% names(DA_MetaData)) DA_MetaData$DiseaseStage = DA_MetaData$DiseaseState
}

Meta_complete = cbind(DA_MetaData, MetaData_extra)
DA_MetaData$SampleID[grepl("^[[:digit:]]+", DA_MetaData$SampleID)] = str_c("S", DA_MetaData$SampleID[grepl("^[[:digit:]]+", DA_MetaData$SampleID)])
write.csv(DA_MetaData, file.path(output_path, "Sample_Info.csv"), row.names=F)

Meta_complete$SampleID[grepl("^[[:digit:]]+", Meta_complete$SampleID)] = str_c("S", Meta_complete$SampleID[grepl("^[[:digit:]]+", Meta_complete$SampleID)])
write.csv(Meta_complete, file.path(output_path, "Sample_Info_complete.csv"), row.names=F)

## Expression & Counts ------

estCount = readRDS(file.path(config$output, paste0(config$prj_name, "_estCount.rds")))

Expression = 2^data_wide - config$count_prior

Sample_list = names(estCount)
Sample_list[grepl("^[[:digit:]]+", Sample_list)] = str_c("S", Sample_list)[grepl("^[[:digit:]]+", Sample_list)]

Exp = data.frame(Gene = rownames(Expression), Expression)
names(Exp) = c("Gene", Sample_list)
Count = data.frame(Gene = rownames(estCount), estCount)
names(Count) = names(Exp)
fwrite(Exp, file.path(output_path, "Gene_Expression_Data.csv"))
fwrite(Count, file.path(output_path, "Gene_Count.csv"))

## comparison info ------
ComparisonID=rownames(comp_info)
c_info=data.frame(ComparisonID, ProjectName=project_info$ProjectID[1], ComparisonContrast = str_c(comp_info$Group_test, ".vs.", comp_info$Group_ctrl), Case_SampleIDs='', Control_SampleIDs='')

unmapped_grp_name_list = c("Phenotype", "population", "SubjectTreatment")

unmapped_name_index = 1
unmapped_record = data.frame("ori"= character(), "new" = character())

for (i in 1:nrow(c_info)) {
  S_meta_sub = MetaData
  Subset_group = comp_info$Subsetting_group[i]
  if (!is.na(Subset_group) && !Subset_group == "") {
    Subset = subset_data(Subset_group, MetaData, estCount)
    S_meta_sub = Subset$S_meta
  }
  DA_meta_sub = DA_MetaData[rownames(S_meta_sub),]

  if (comp_info$Analysis_method[i] == "DESeq2") {
    c_info$ComparisonID[i] = paste0(c_info$ComparisonID[i], "_DESeq")
  } else if (comp_info$Analysis_method[i] == "limma") {
    c_info$ComparisonID[i] = paste0(c_info$ComparisonID[i], "_limma")
  }    
  compare_group_ori = comp_info$Group_name[i]
  print(compare_group_ori)
  print(head(S_meta_sub))
  c_info$Case_SampleIDs[i]=str_c(DA_meta_sub[S_meta_sub[, compare_group_ori]==comp_info$Group_test[i], 'SampleID'], collapse=";") 
  c_info$Control_SampleIDs[i]=str_c(DA_meta_sub[S_meta_sub[, compare_group_ori]==comp_info$Group_ctrl[i], 'SampleID'], collapse=";") 
  
  compare_group = header_mapping$DA_header[header_mapping$NGSone_header == compare_group_ori]

  if (length(compare_group) > 0 && !compare_group %in% c("SampleType", 'Tissue', "Ethnicity")) {
    c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
    c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
  } else if (length(compare_group) > 0 && compare_group %in% c('Tissue', "Ethnicity")) {
    if (comp_info$Group_test[i] %in% DA_meta_sub[, compare_group] && comp_info$Group_ctrl[i] %in% DA_meta_sub[, compare_group]) {
      c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
      c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
    } else if (unmapped_name_index < 4 && !compare_group_ori %in% unmapped_record$ori) {
      compare_group = unmapped_grp_name_list[unmapped_name_index]
      unmapped_record[unmapped_name_index,] = c(compare_group_ori, compare_group)
      unmapped_name_index = unmapped_name_index + 1
      c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
      c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
    } else if (unmapped_name_index < 4 && compare_group_ori %in% unmapped_record$ori) {
      compare_group = unmapped_record$new[which(unmapped_record$ori == compare_group_ori)]
      c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
      c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
    }
  } else if (length(compare_group) == 0 && unmapped_name_index < 4) {
    if (!compare_group_ori %in% unmapped_record$ori) {
      compare_group = unmapped_grp_name_list[unmapped_name_index]
      unmapped_record[unmapped_name_index,] = c(compare_group_ori, compare_group)
      unmapped_name_index = unmapped_name_index + 1
      c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
      c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
    } else if (compare_group_ori %in% unmapped_record$ori) {
      compare_group = unmapped_record$new[which(unmapped_record$ori == compare_group_ori)]
      c_info[i, paste0("Case_", compare_group)] = comp_info$Group_test[i]
      c_info[i, paste0("Control_", compare_group)] = comp_info$Group_ctrl[i]
    }
  }

  compare_group_new = compare_group
  if (compare_group == "Age") {
    compare_group_new = "AgeCategory"
  } else if (compare_group == "SampleSource" && "Cell_Line" %in% meta_names_ori) {
    compare_group_new = "CellLine"
  } else if (compare_group == "Transfection" && "Genotype" %in% meta_names_ori) {
    compare_group_new = "Genotype"
  } else if (compare_group == "SampleType" && "group" %in% meta_names_ori) {
    compare_group_new = "SubjectGroup"
  }
  if (! compare_group_new == compare_group) {
    c_info[i, paste0("Case_", compare_group_new)] = comp_info$Group_test[i]
    c_info[i, paste0("Control_", compare_group_new)] = comp_info$Group_ctrl[i]
  }

  Covariate_list = setdiff(names(DA_meta_sub), c("SampleID", "ProjectName", "Organism", compare_group))
  
  if (length(Covariate_list) > 0) {
    for (j in 1:length(Covariate_list)) {
      case = unique(DA_meta_sub[S_meta_sub[, compare_group_ori]==comp_info$Group_test[i], Covariate_list[j]])
      control = unique(DA_meta_sub[S_meta_sub[, compare_group_ori]==comp_info$Group_ctrl[i], Covariate_list[j]])

      Covariate = Covariate_list[j]
      if (Covariate_list[j] == "Age") {
        Covariate = "AgeCategory"
      } else if (Covariate_list[j] == "SampleSource" && "Cell_Line" %in% meta_names_ori) {
        Covariate = "CellLine"
      } else if (Covariate_list[j] == "Transfection" && "Genotype" %in% meta_names_ori) {
        Covariate = "Genotype"
      } else if (Covariate_list[j] == "SampleType" && "group" %in% meta_names_ori) {
        Covariate = "SubjectGroup"
      }
      
      if (length(case)== 1) {
        c_info[i, paste0("Case_", Covariate)]=case
        if (Covariate_list[j] %in% c("SampleSource", 'Transfection')) {
          c_info[i, paste0("Case_", Covariate_list[j])]=case
        }
      }
      if (length(control)== 1) {
        c_info[i, paste0("Control_", Covariate)]=control
        if (Covariate_list[j]  %in% c("SampleSource", 'Transfection')) {
          c_info[i, paste0("Control_", Covariate_list[j])]=control
        }
      }
    }
  }
} 

if (nrow(unmapped_record) > 0) {
  for (i in 1:nrow(c_info)) {
    for (j in 1: nrow(unmapped_record)) {
      case = unique(S_meta_sub[S_meta_sub[, comp_info$Group_name[i]]==comp_info$Group_test[i], unmapped_record$ori[j]])
      control = unique(S_meta_sub[S_meta_sub[, comp_info$Group_name[i]]==comp_info$Group_ctrl[i], unmapped_record$ori[j]])
      if (length(case)== 1) c_info[i, paste0("Case_", unmapped_record$new[j])]=case
      if (length(control)== 1) c_info[i, paste0("Control_", unmapped_record$new[j])]=control
    }
  }
}

fwrite(c_info, file.path(output_path, "Comparisons_Info.csv"))


## comparison data ------
comp_data = data_results[, which(str_detect(names(data_results), 
								"log2FoldChange|pvalue|padj|logFC|P.value|Adj.P.value"))]
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
