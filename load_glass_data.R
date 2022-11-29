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


readcounts <- readcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% metadata$GS_ID) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()


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


cell_cycle <- read_csv("cell_cycle.csv")
 
 DNA_binding <- read_csv("DNA_binding.csv")
 
 primrecurvolcano <- deseq2Resultsprimrecur %>% 
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
 dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colors)) %>% 
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>% 
   dplyr::mutate(size = ifelse(colors == "DNA binding" | colors == "Cell cycle", 1.5,1))
 
 ggplot(primrecurvolcano, aes(log2FoldChange, -log10(padj))) + 
   geom_point(aes(col=colors, size = size)) + 
   scale_radius(range = c(1,2))+ 
   guides(size = "none")+
   scale_color_manual(values=c("purple","royalblue", "gray87", "gray40")) + 
   theme_bw()+
   geom_text_repel(data=subset(primrecurvolcano, colors == "Cell cycle" | colors == "DNA binding"), aes(label=gene_name), max.overlaps = 4) +
   geom_hline(yintercept=1.3, linetype="dashed")+
   geom_vline(xintercept = 0.75, linetype="dashed")+
   geom_vline(xintercept = -0.75, linetype="dashed")+
   labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between primary (n = 71) and recurrent (n= 69) samples')+
   xlab('Log fold change')+
   ylab('-log10(adjusted p-value)')+
   xlim(-7,7)
# 
# 
# 
signupRNAprimrecur <- deseq2Resultsprimrecur %>% 
   dplyr::filter(padj <0.05) %>% 
   dplyr::filter(log2FoldChange > 0.75) %>% 
   left_join(gene.annot2)
 
 
 #948 genes up 
signdownRNAprimrecur <- deseq2Resultsprimrecur %>% 
   dplyr::filter(padj <0.05) %>% 
   dplyr::filter(log2FoldChange < (-0.75)) %>% 
   left_join(gene.annot2)


signdownRNAprimrecurCOMP <- signdownRNAprimrecur %>% 
  inner_join(signdownRNA, by = "gene_name")


143/173

signupRNAprimrecurCOMP <- signupRNAprimrecur %>% 
  inner_join(signupRNA, by = "gene_name")

790/948

`%!in%` = Negate(`%in%`)

metadata <- metadata %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Sample_Name))

metadata$Sample_Number[!duplicated(metadata$Sample_Number)]
  
duplicates <- metadata$Sample_Number[duplicated(metadata$Sample_Number)]

sample_unique <- metadata %>% 
  dplyr::filter(Sample_Number %!in% duplicates)

metadata_new <- metadata %>% 
  dplyr::filter(Sample_Type == "R")


sum(duplicated(metadata_new$GLASS_ID))

metadata_new <- metadata_new[duplicated(metadata_new$GLASS_ID),]

metadata_onlylastresec <- metadata %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_005" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_006" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_007" & resection == "S3")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_007" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_029" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_032" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_032" & resection == "S3")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_043" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_045" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_133" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_147" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_163" & resection == "S3")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_171" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_172" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_AUMC_034" & resection == "S2")) %>% 
  dplyr::filter(!(GLASS_ID == "GLNL_EMCR_129" & resection == "S2"))

metadata_onlylastresec[duplicated(metadata_onlylastresec$GLASS_ID),] 

metadata_onlylastresecRec <-  metadata_onlylastresec%>% dplyr::filter(Sample_Type == "R")

readcounts_lastresec <- readcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% metadata_onlylastresec$GS_ID) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()

dds_lastres <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts_lastresec,
                                      colData = metadata_onlylastresec,
                                      design= ~Sample_Type)
dds_lastres <-DESeq(dds_lastres)
vsd_lastres <- DESeq2::vst(dds_lastres)
readcount.vst_lastres <- assay(vsd) %>% as.data.frame 



deseq2Resultsprimrecur_lastres <- results(dds_lastres) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)

rm(vsd,tmp, dds)

#867 genes up
signupRNAprimrecur_lastres <- deseq2Resultsprimrecur_lastres %>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(log2FoldChange > 0.75) 

#209 genes down 
signdownRNAprimrecur_lastres <- deseq2Resultsprimrecur_lastres %>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(log2FoldChange < (-0.75))


signupRNAprimrecurCOMP_lastres <- signupRNAprimrecur_lastres %>% 
  inner_join(signupRNA, by = "gene_name")

signdownRNAprimrecurCOMP_lastres <- signdownRNAprimrecur_lastres %>% 
  inner_join(signdownRNA, by = "gene_name")





primrecurvolcano_lastres <- deseq2Resultsprimrecur_lastres %>% 
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% HOXgenes$gene_name, "HOX gene", colors)) %>% 
  dplyr::mutate(size = ifelse(colors == "HOX gene", 3,1))

#%>% 
#  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colors)) %>% 
#  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>% 
#  dplyr::mutate(size = ifelse(colors == "DNA binding" | colors == "Cell cycle", 1.5,1))


pdf("Volcano_RNA_primrecur.pdf", width = 10, height = 6)
ggplot(primrecurvolcano_lastres, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col=colors, size = size)) + 
  scale_radius(range = c(1,2))+ 
  geom_text_repel(data=subset(primrecurvolcano_lastres, colors == "Significant" | colors == "HOX gene"), aes(label=gene_name), max.overlaps = 7) +
  guides(size = "none")+
  scale_color_manual(values=c("royalblue", "gray87", "gray40")) + 
  theme_bw()+
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between primary (n = 71) and recurrent (n= 79) samples')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-7,7)



dev.off()

