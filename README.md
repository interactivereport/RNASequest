# RNASequest: An end-to-end reproducible RNAseq data analysis and publishing frame-work

![RNASequest](https://interactivereport.github.io/RNASequest/Figure1_sm.png?raw=true "RNASequest")

**Fig. 1.** Overview of the RNASequest workflow. (A) Analyst supplies gene expression matrix, sample information to the automated pipeline, ExpressionAnalysis, in abbreviation EA, to generate reports in Bookdown and interactive slide deck formats, and a data visualization app for Biologists with limited computational experience to investigate datasets by reviewing the reports and exploring the data interactively. (B) EA will check covariates and guide ana-lysts to build correct models for differential gene expression analysis. (C) R data objects outputted by EA will be uploaded to Quickomics R Shiny application for further exploration and visualization in PCA, Heatmap, Pathway, Volcano, Boxplot, and Venn Diagram. (D) EA publish module will automatically generate analysis report in both Bookdown and interactive online slides format.  (E) ShinyOne, a R Shiny app will manage the collection of datasets with Quickomics launching links and links to Bookdown documents and slide decks. It provides basic search and sorting functions for users to locate datasets of interest. 

## Expression analysis (EA) component

A pipeline to analysis RNAseq

Four main functions are provided:

  - EAinit: Generate a set of project analysis files based on a DNAnexus result folder;
  - EAqc: Analyze the covariates against the expression to determine if the expression is needed to be adjusted;
  - EArun: Produce QuickOmics object for webserver loading.
  - EA2DA: Produce required data files for DiseaseAtlas project import.

### Installation/Set up
Add "/camhpc/ngs/tools/expressionAnalysis" into your **_PATH_** environment variable.

### EAinit
```
EAinit A/path/to/a/DNAnexus/result/folder
```
The execution of above command will create a sub-folder (QuickOmics_[timestamp]) in the specified DNAnexus result folder.
There will be five files in the folder:

- compareInfo.csv: an empty comparison definition file (with header). Please fill in this file before ```EArun``` call.
- config.yml: a config fill specifies the parameters of the ```EAqc``` and ```EArun```. Please update **covariates_adjust** after ```EAqc```.
- geneAnnotation.csv: a gene annotation file including gene symbol.
- sampleMeta.csv: a sample meta information file, please feel free to add additional columns whose column names should be considered to be added into **covariates_check** in *config.yml*.
- alignQC.pdf: plots generated from alignment QC metrics.

**_Please pay attention on the std out messages._**

### EAqc
```
EAqc A/path/to/a/config/file
```
The execution of the command with the above default config file, expression PC analysis will be done against covariates specified in **covariates_check** in *config.yml* file. An excel file will list p-value for all numeric and categorical covariates, and with significant ones will be in plot pdf files. The analysis before covariate adjusting will have prefix *covariatePCanalysis_noAdjust*. 
Based on the above results, you can add covariates into **covariates_adjust** in *config.yml* file, and run ```EAqc``` again. This time additional expression PC analysis will be applied to covariate adjusted expression with files started with *covariatePCanalysis_Adjusted*. 

**_Please pay attention on the std out messages._**

### EArun
```
EArun A/path/to/a/config/file
```
**_Please fill the compareInfo.csv before executing the above command_**

The execution of above command will produce R object for QuickOmics webserver to load. The process will generate the covariate adjusted logTPM for visualization; complete differentially expressed gene analysis and gene network generation. 

The results (three files) currently will need to be copied to a folder on ngs, in order to have web accession. **_Please pay attention on the std out messages._**

### EA2DA
```
EA2DA A/path/to/a/config/file
```
The execution of above command will produce 6 data files which are required for the DiseaseAtlas project import.

**_Please fill the empty entries in the Project_Info.csv before import**

### Administration
There are two config files in the pipeline folder:
 - config.tmp.yml: The template of the config file, with all default values;
 - sys.yml: the system config file, which includes:
    1. genome_path: the root path where the genome definition files (gtf) are located
    2. notCovariates: the column names from the sample meta information should not be considred as default covariates
    3. qc2meta: the column names from mapping QC file should be extracted and inserted into sample meta table
    4. QuickOmics_path: the file path to store the files for QuickOmics web server display
    4. DA_columns: the column names available for the sample meta table in the DiseaseAtlas system

## Quickomic component

GitHub: https://github.com/interactivereport/Quickomics

Tutorial: https://interactivereport.github.io/Quickomics/tutorial/docs/

## Bookdown component

https://interactivereport.github.io/RNASequest/tutorial/docs/bookdown-report.html#bookdown-report

## Slide deck component

https://interactivereport.github.io/RNASequest/tutorial/docs/online-slide-deck.html#online-slide-deck


## ShinyOne component

https://interactivereport.github.io/RNASequest/tutorial/docs/shinyone.html#shinyone


