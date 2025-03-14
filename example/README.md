# An example using "recount3"

The "recount3.example.R" is an example code of using the "recount3" package to retrieve one public dataset (GSE131869) by SRA record (SRP199678) to be used by *RNAsequest*. If it was executed, the folder "SRP199678" will be created, along with count matrix, sample meta data, gene annotation file, et al. to be used by **EAinit** from *RNAsequest*.

After executing **EAinit** with the path to the above "SRP199678" folder, an RNAsequest folder "EA........_X" x will be created where the project analysis will be housed. "........" is the date as YYYYMMDD, and the "X" is a numeric suffix, e.g. "EA20230103_0". The "config.yml" contains the project information which is the sole parameter to **EAqc** and **EArun**, and the "data" folder contains input data used by *RNAsequest* pipeline, including count matrix, sample meta data and comparison definition file. 
