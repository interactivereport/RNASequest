# Expression analysis (EA)

A pipeline to analysis RNAseq

Four main functions are provided:

  - EAinit: Generate a set of project analysis files based on a DNAnexus result folder;
  - EAqc: Analyze the covariates against the expression to determine if the expression is needed to be adjusted;
  - EArun: Produce QuickOmics object for webserver loading.
  - EA2DA: Produce required data files for DiseaseAtlas project import.

# Installation/Set up
Add "/camhpc/ngs/tools/expressionAnalysis" into your **_PATH_** environment variable.

# EAinit
```
EAinit A/path/to/a/DNAnexus/result/folder
```
The execution of above command will create a sub-folder (QuickOmics_[timestamp]) in the specified DNAnexus result folder.
There will be four files in the folder:
    - compareInfo.csv: an empty comparison definition file (with header). Please fill in this file before ```EArun``` call.
    - config.yml: a config fill specifies the parameters of the ```EAqc``` and ```EArun```. Please update **covariates_adjust** after ```EAqc```.
    - geneAnnotation.csv: a gene annotation file including gene symbol.
    - sampleMeta.csv: a sample meta information file, please feel free to add additional columns whose column names should be considered to be added into **covariates_check** in *config.yml*.

**_Please pay attention on the std out messages._**

# EAqc
```
EAqc A/path/to/a/config/file
```
The execution of the command with the above default config file, expression PC analysis will be done against covariates specified in **covariates_check** in *config.yml* file. An excel file will list p-value for all numeric and categorical covariates, and with significant ones will be in plot pdf files. The analysis before covariate adjusting will have prefix *covariatePCanalysis_beforeCovariateAdjust*. 
Based on the above results, you can add covariates into **covariates_adjust** in *config.yml* file, and run ```EAqc``` again. This time additional expression PC analysis will be applied to covariate adjusted expression with files started with *covariatePCanalysis_afterCovariateAdjust*. 

**_Please pay attention on the std out messages._**

# EArun
```
EArun A/path/to/a/config/file
```
**_Please fill the compareInfo.csv before executing the above command_**

The execution of above command will produce R object for QuickOmics webserver to load. The process will generate the covariate adjusted logTPM for visualization; complete differentially expressed gene analysis and gene network generation. 

The results (three files) currently will need to be copied to a folder on ngs, in order to have web accession. **_Please pay attention on the std out messages._**

# EA2DA
```
EA2DA A/path/to/a/config/file
```
The execution of above command will produce 6 data files which are required for the DiseaseAtlas project import.

**_Please fill the empty entries in the Project_Info.csv before import**

# Administration
There are two config files in the pipeline folder:
 - config.tmp.yml: The template of the config file, with all default values;
 - sys.yml: the system config file, which includes:
    1. genome_path: the root path where the genome definition files (gtf) are located
    2. notCovariates: the column names from the sample meta information should not be considred as default covariates
    3. qc2meta: the column names from mapping QC file should be extracted and inserted into sample meta table
    4. QuickOmics_path: the file path to store the files for QuickOmics web server display
    4. DA_columns: the column names available for the sample meta table in the DiseaseAtlas system





