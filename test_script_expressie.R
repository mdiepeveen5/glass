library(tidyverse)
library(DESeq2)
library(ggplot2)

tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))

# patient met read counts < 750.000 eruit scoopt en opnieuw

# DESeq2 - volg manual en maak PCA plot van Glass RNA data

# load RNAseq ----

readcounts <- tmp %>% 
  dplyr::mutate(Chr=NULL) %>% 
  dplyr::mutate(Start=NULL) %>% 
  dplyr::mutate(Length=NULL) %>% 
  dplyr::mutate(End=NULL) %>% 
  dplyr::mutate(Strand=NULL) %>% 
  tibble::column_to_rownames("Geneid") %>% 
  dplyr::mutate(sumrow = rowSums(.)) %>% 
  dplyr::filter(sumrow > 230 * 3)%>% 
  dplyr::mutate(sumrow=NULL) %>% 
  t() %>% 
  as.data.frame() %>%  
  dplyr::mutate(rs = rowSums(.)) %>% 
  tibble::rownames_to_column('GS_ID') %>% 
  dplyr::mutate("GS_ID" = gsub(".","-", GS_ID, fixed = T)) %>% 
  dplyr::filter(rs >750000) %>% 
  dplyr::mutate(rs =NULL) %>% 
  tibble::column_to_rownames('GS_ID') %>% 
  t() %>% 
  as.data.frame()
  
# metadata ----

metadata<-
  data.frame(GS_ID = colnames(readcounts)) %>% 
  dplyr::left_join(
    read.csv("~/projects/glass/data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") ,
    by=c('GS_ID'='GS_ID')
  ) %>% 
  # dplyr::mutate(Sample_Type=NULL) %>% 
  # dplyr::mutate(Recurrent_Type=NULL) %>% 
  dplyr::mutate(Exclude=NULL) %>%   
  dplyr::mutate(Surgery_ID=NULL) %>% 
  dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex))

# supervised analysis ----

## deseq2 ----

#stopifnot(colnames(totaldata) == coldata$GS_ID)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts,
                                      colData = metadata,
                                      design= ~)
dds <-DESeq(dds)
vsd <- DESeq2::vst(dds, blind=T)

res <-results(dds) %>% 
  as.data.frame %>% 
  dplyr::arrange(padj,pvalue)

selected.features.supervised <- res %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::select('gene_id') %>% 
  dplyr::slice_head(n=1000)




readcount.vst <- assay(vsd) %>% as.data.frame 

selected.features <- readcount.vst %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=1000) %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::select('gene_id')

readcount.vst <-readcount.vst %>% 
  tibble::rownames_to_column('gene_id')
  
total <- dplyr::left_join(selected.features, readcount.vst, by = "gene_id") %>% 
  tibble::column_to_rownames('gene_id')
  

supervised <- dplyr::left_join(selected.features.supervised, readcount.vst, by = "gene_id") %>% 
  tibble::column_to_rownames('gene_id')

#gen namen sorteden op padj
pca <- total %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')

pca_supervised <- supervised %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')
  

# combi <- total %>%
#  dplyr::left_join(pca$x) %>%
#  dplyr::filter(padj < 0.01)

#plotPCA(vsd, intgroup = c("Sample_Sex"))

unsupervised_sorted <-dplyr::left_join(metadata, pca)
supervised_sorted <-dplyr::left_join(metadata, pca_supervised)

ggplot(unsupervised_sorted, aes(x=PC3, y=PC4, color = Sample_Sex)) + geom_point()

ggplot(supervised_sorted, aes(x=PC1, y=PC2, color = Sample_Sex, group = GLASS_ID)) +
  geom_point() +
  geom_line() +
  theme_classic()


