suppressWarnings(suppressMessages(require(tidyverse)))
suppressMessages(require(data.table))
suppressMessages(require(GenomicRanges))
suppressMessages(require(R.utils))

gtf.gz_to_gene_info <- function(strGTF,strOut=NULL){
    #I have tested the script with the following files
    #Enter gtf file name
    #file="/camhpc/ngs/projects/DNAnexus_references/rnaseq/human/Human.GRCh38.v34.l1_5.ERCC/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf.gz"
    #file="/camhpc/ngs/projects/DNAnexus_references/rnaseq/mouse/Mouse.GRCm38.102/Mus_musculus.GRCm38.transcript.gtf.gz"
    #file="/camhpc/ngs/projects/DNAnexus_references/rnaseq/monkey/Monkey.Macaca_fascicularis_5.0_NCBI_RefSeq.aavhSMN1.REV4_2/Monkey_REV4b.transcript.gtf.gz"
    #file="/camhpc/ngs/projects/DNAnexus_references/rnaseq/monkey/Monkey.macFas5.v100.l1_5/Monkey.macFas5.v100.l1_5.transcript.transcript.gtf.gz"
     #file="/camhpc/ngs/projects/DNAnexus_references/rnaseq/rat/Rattus_norvegicus_ensemble_Rnor6.0.89/Rattus_norvegicus_ensemble_Rnor6.0.89.transcript.transcript.gtf.gz"

    gene_info_file <- strOut
    #if(is.null(gene_info_file)) gene_info_file=str_replace(strGTF, "gtf.gz", "gene_info.csv")  #will save gene information file here
    if(is.null(gene_info_file)) gene_info_file=str_replace(str_replace(strGTF, "gtf.gz", "gene_info.csv"),"gtf$","gene_info.csv") #O'Young
    if(gene_info_file==strGTF) gene_info_file <- paste0(strGTF,".gene_info.csv")#O'Young
    #gtf<-fread(strGTF)
    if ( str_detect(strGTF, "gtf.gz$") ) {  
        gtf=fread(cmd=str_c('zcat ', strGTF, '| grep -v "^#" ') ) #use zcat for gtf.gz files
    } else {
        gtf<-fread(cmd=paste0('grep -v "^#" ',strGTF)) #O'Young
    }
    gtf%>%group_by(V3)%>% (dplyr::count)
    #some files use gene_biotype, some files use gene_type
    if (sum(str_detect(gtf$V9, "gene_biotype"))>100) {gene_type_key="gene_biotype"} else {gene_type_key="gene_type"}
    
    if (sum(gtf$V3=="transcript")>100) { #standard annotations have transcript
        tx_info<-gtf%>%filter(V3=="transcript")%>%tidyr::extract(V9, c("geneID"), 'gene_id "(.+?)";', remove=F) %>%
            tidyr::extract(V9, c("gene_name"), 'gene_name "(.+?)";', remove=F)%>%
            tidyr::extract(V9, c("gene_type"), str_c(gene_type_key, ' "(.+?)";'), remove=F)
    } else { #some annotation files only have exons, like Monkey_REV4b.transcript.gtf.gz
        tx_info<-gtf%>%filter(V3=="exon")%>%tidyr::extract(V9, c("geneID"), 'gene_id "(.+?)";', remove=F) %>%
            tidyr::extract(V9, c("gene_name"), 'gene_name "(.+?)";', remove=F)%>%
            tidyr::extract(V9, c("gene_type"), str_c(gene_type_key, ' "(.+?)";'), remove=F)
    }
    
    colSums(is.na(tx_info)) #if there are NAs, check why the information is not extracted. Some files don't have gene type info.
    
    
    
    #add ERCC. If no ERCC, get empty data.table which won't affect the results
    ERCC<-gtf%>%filter(str_detect(V1, "ERCC"))%>%tidyr::extract(V9, c("geneID"), 'gene_id "(.+?)";', remove=F) %>%
        tidyr::extract(V9, c("gene_name"), 'gene_name "(.+?)";', remove=F)%>%
        tidyr::extract(V9, c("gene_type"), str_c(gene_type_key, ' "(.+?)";'), remove=F)
    #add plasmid genes. If no plasmid genes, get empty data.table which won't affect the results
    plasmid_genes<-gtf%>%filter(V2=="gb2gtf")%>%tidyr::extract(V9, c("geneID"), 'gene_id "(.+?)";', remove=F) %>%
        tidyr::extract(V9, c("gene_name"), 'transcript_id "(.+?)";', remove=F)%>%
        tidyr::extract(V9, c("gene_type"), 'Feature_type "(.+?)"', remove=F)
    
    all_genes<-rbind(ERCC%>%filter(!duplicated(geneID))%>%select(geneID, gene_name, gene_type),
    plasmid_genes%>%filter(!duplicated(geneID))%>%select(geneID, gene_name, gene_type),
    tx_info%>%filter(!duplicated(geneID))%>%select(geneID, gene_name, gene_type) ) %>%filter(!duplicated(geneID))
    
    #now get gene length by combining all exons for a gene
    exons<-gtf%>%filter(V3=="exon")%>%tidyr::extract(V9, c("geneID"), 'gene_id "(.+?)";', remove=F)
    gr1=GRanges(seqnames=exons$V1,  ranges=IRanges(exons$V4, exons$V5), strand=exons$V7, geneID=exons$geneID)
    gr2=split(gr1, exons$geneID)  #group exons by geneID
    gr2=GenomicRanges::reduce(gr2) #merge overlapping exons
    gr2=unlist(gr2, use.names=T)
    info1=data.frame(geneID=names(gr2), width=width(gr2) )
    Gene_Length<-info1%>%group_by(geneID)%>%summarize(Length=sum(width))
    all_genes<-all_genes%>%left_join(Gene_Length)
    colSums(is.na(all_genes)); dim(all_genes) #check to make sure thing look all right
    fwrite(all_genes, gene_info_file)
}
