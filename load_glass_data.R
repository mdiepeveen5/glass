library(tidyverse)
library(DESeq2)
library(ggplot2)
#---- load data

tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))

# load RNAseq ----
# patient met read counts < 750.000 eruit scoopt en opnieuw

readcounts <- tmp %>% 
  dplyr::mutate(Chr=NULL) %>% 
  dplyr::mutate(Start=NULL) %>% 
  dplyr::mutate(Length=NULL) %>% 
  dplyr::mutate(End=NULL) %>% 
  dplyr::mutate(Strand=NULL) %>% 
  dplyr::filter(Geneid %in% c("ENSG00000280987.4", "ENSG00000187522.16", "ENSG00000285053.1")== F) %>% 
  tibble::column_to_rownames("Geneid") %>% 
  dplyr::mutate(sumrow = rowSums(.)) %>% 
  dplyr::filter(sumrow > 229 * 3)%>% 
  dplyr::mutate(sumrow=NULL) %>% 
  t() %>% 
  as.data.frame() %>%  
  dplyr::mutate(rs = rowSums(.)) %>% 
  tibble::rownames_to_column('GS_ID') %>% 
  dplyr::mutate("GS_ID" = gsub(".","-", GS_ID, fixed = T)) %>% 
  dplyr::filter(rs >750000) %>% 
  dplyr::mutate(rs =NULL) %>% 
  dplyr::filter(GS_ID %in% metadata$GS_ID) %>% 
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
  dplyr::mutate(Exclude=NULL) %>%   
  dplyr::mutate(Surgery_ID=NULL) %>% 
  dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex)) %>% 
  dplyr::left_join(metadataALL) %>% 
  dplyr::filter(is.na(Sample_Type) == F) %>% 
  dplyr::filter(Sample_Type != "X")


#readcountsPR1 <- readcounts %>% 
#  t() %>% 
#  as.data.frame() %>% 
#  rownames_to_column("GS_ID") %>% 
#  dplyr::filter(GS_ID %in% metadata$GS_ID) %>% 
#  column_to_rownames("GS_ID") %>% 
#  t() %>% 
#  as.data.frame()


# vst transformed counts ----

dds <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts,
                                      colData = metadata,
                                      design= ~Sample_Type)
dds <-DESeq(dds)
vsd <- DESeq2::vst(dds)
readcount.vst <- assay(vsd) %>% as.data.frame 



deseq2Resultsprimrecur <- results(dds) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)

rm(vsd,tmp, dds)


primrecurvolcano <- deseq2Resultsprimrecur %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colours)) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colours)) %>% 
  dplyr::mutate(size = ifelse(colours == "DNA binding" | colours == "Cell cycle", 1.5,1))

ggplot(primrecurvolcano, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col=colours, size = size)) + 
  scale_radius(range = c(1,2))+ 
  guides(size = "none")+
  scale_color_manual(values=c("purple","royalblue", "gray87", "gray40")) + 
  theme_bw()+
  geom_text_repel(data=subset(primrecurvolcano, colours == "Cell cycle" | colours == "DNA binding"), aes(label=gene_name), max.overlaps = 4) +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between primary (n = 71) and recurrent (n= 69) samples')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-7,7)



deseq2ResultsSIGNpr <- primrecurvolcano%>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  left_join(gene.annot2)


#173 genes down
signupRNAprimrecur <- deseq2ResultsSIGNpr %>% 
  dplyr::filter(log2FoldChange >0)

#948 genes up 
signdownRNAprimrecur <- deseq2ResultsSIGNpr %>% 
  dplyr::filter(log2FoldChange <0)
