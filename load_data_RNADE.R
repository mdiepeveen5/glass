library(tidyverse)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(readr)
library(readODS)

#---- load data
cycling.cell.markers <- readODS::read_ods('scRNA-cycling-genes.ods',col_names=F) %>%
  dplyr::mutate(A = ifelse(duplicated(A),NA,A)) %>%
  dplyr::mutate(B = ifelse(duplicated(B),NA,B)) %>%
  dplyr::mutate(C = ifelse(duplicated(C),NA,C)) %>%
  dplyr::mutate(E = ifelse(duplicated(E),NA,E)) %>%
  dplyr::mutate(F = ifelse(duplicated(F),NA,F)) %>%
  dplyr::rename(G1.S.tirosh = A) %>%
  dplyr::rename(G2.M.tirosh = B) %>%
  dplyr::rename(cycling.melanoma.tirosh = C) %>%
  dplyr::rename(G1.S.neftel = E) %>%
  dplyr::rename(G2.M.neftel = F) %>%
  dplyr::mutate(D = NULL) %>%
  dplyr::mutate(G = NULL) %>%
  dplyr::mutate(H = NULL) %>%
  tidyr::pivot_longer(cols=colnames(.)) %>%
  tidyr::drop_na(value) %>%
  dplyr::filter(value %in% c('10.1126/science.aad0501','10.1016/j.cell.2019.06.024','G1/S phase-specific','G2/M phase-specific','melanoma cell cycle genes','G1/S phase-specific','G2/M phase-specific') == F) %>%
  dplyr::mutate(marker=T) %>%
  dplyr::mutate(name = as.factor(name)) %>%
  dplyr::mutate(value = as.factor(value)) %>%
  dplyr::rename(source = name) %>%
  dplyr::rename(gene_name = value) %>%
  pivot_wider(names_from = source, values_from = marker, values_fill = F) %>%
  dplyr::mutate(gene_name = as.character(gene_name))

load(file = "metadataALL.Rdata")

gene.annot <- read.table('data/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf', comment.char='#',sep="\t",header=F)
gene.annot2 <- gene.annot %>% 
  dplyr::filter(V3 =='gene') %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id ([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name([^;]+);.*$","\\1",V9)) %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::mutate(gene_name = gsub("^[ ]+","",gene_name))

rm(gene.annot)
tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))

# load RNAseq ----
#alleen patienten met A_IDH en IDH_HG plus genen met readcounts > 230*3
#firsts run readcountsDE without GS_ID %in% metadataDE, then run metadataDE, then run readcountsDE with GS_ID %in% metadataDE
readcountsDE <- tmp %>% 
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


metadataDE<-
  data.frame(GS_ID = colnames(readcountsDE)) %>% 
  dplyr::left_join(
    read.csv("~/projects/glass/data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") ,
    by=c('GS_ID'='GS_ID')
  ) %>% 
  dplyr::mutate(Exclude=NULL) %>%   
  dplyr::mutate(Surgery_ID=NULL) %>% 
  dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex)) %>% 
  dplyr::left_join(metadataALL) %>% 
  dplyr::filter(methylation.sub.diagnosis %in% c("A_IDH_HG", "A_IDH")) 

readcountsDE <- readcountsDE %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% metadataDE$GS_ID) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()

saveRDS(metadataDE, file = "/tmp/proteomics/metadataDE.Rds")
saveRDS(readcountsDE, file = "/tmp/proteomics/readcountsDE.Rds")
  

stopifnot(length(colnames(readcountsDE)) == length(metadataDE$GS_ID))

# vst transformed counts ----

dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE,
                                      colData = metadataDE,
                                      design= ~methylation.sub.diagnosis)
dds2 <-DESeq(dds2)
vsd2 <- DESeq2::vst(dds2)
readcount.vstDE <- assay(vsd2) %>% as.data.frame 

save(readcount.vstDE, file = "readcount.vstDE.Rdata")

deseq2Results <- results(dds2) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)


save(deseq2Results, file = "deseq2Results.Rdata")


deseq2ResultsSIGN <- deseq2Results%>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  left_join(gene.annot2)

signupRNA <- deseq2ResultsSIGN %>% 
  dplyr::filter(log2FoldChange >0)

signdownRNA <- deseq2ResultsSIGN %>% 
  dplyr::filter(log2FoldChange <0)

upprotRNA <- signupRNA %>% 
  inner_join(diffprot2signup, by = c("gene_name"= "Protein"))

downprotRNA <- signdownRNA %>% 
  inner_join(diffprot2signdown, by = c("gene_name"= "Protein"))


summary(deseq2Results)
plotMA(deseq2Results)

rm(vsd2,tmp, dds2, readcount.vst2)


#alleen patienten met A_IDH en IDH_HG plus genen met readcounts > 230*3
# 
# WHOclass0305<- WHOclass0305<- read_csv("WHOclassification_03052022.csv", show_col_types = FALSE) %>% 
#   dplyr::filter(Sample_Name %in% c("007_R4","026_R2", "219_R3", "216_R3") ==F) %>% 
#   dplyr::mutate(Grade = ifelse(WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 4", "Grade4", "Grade2or3"))
# 
# 
# 
# ###GRADE
# readcountsDEGRADE <- tmp %>% 
#   dplyr::mutate(Chr=NULL) %>% 
#   dplyr::mutate(Start=NULL) %>% 
#   dplyr::mutate(Length=NULL) %>% 
#   dplyr::mutate(End=NULL) %>% 
#   dplyr::mutate(Strand=NULL) %>% 
#   dplyr::filter(Geneid %in% c("ENSG00000280987.4", "ENSG00000187522.16", "ENSG00000285053.1")== F) %>% 
#   tibble::column_to_rownames("Geneid") %>% 
#   dplyr::mutate(sumrow = rowSums(.)) %>% 
#   dplyr::filter(sumrow > 229 * 3)%>% 
#   dplyr::mutate(sumrow=NULL) %>% 
#   t() %>% 
#   as.data.frame() %>%  
#   dplyr::mutate(rs = rowSums(.)) %>% 
#   tibble::rownames_to_column('GS_ID') %>% 
#   dplyr::mutate("GS_ID" = gsub(".","-", GS_ID, fixed = T)) %>% 
#   dplyr::filter(rs >750000) %>% 
#   dplyr::mutate(rs =NULL) %>% 
#   dplyr::filter(GS_ID %in% metadataDEGRADE$GS_ID) %>% 
#   tibble::column_to_rownames('GS_ID') %>% 
#   t() %>% 
#   as.data.frame()
# 
# metadataDEGRADE<-
#   data.frame(GS_ID = colnames(readcountsDEGRADE)) %>% 
#   dplyr::left_join(
#     read.csv("~/projects/glass/data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") ,
#     by=c('GS_ID'='GS_ID')
#   ) %>% 
#   dplyr::mutate(Exclude=NULL) %>%   
#   dplyr::mutate(Surgery_ID=NULL) %>% 
#   dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex)) %>% 
#   dplyr::left_join(metadataALL) %>% 
#   dplyr::inner_join(WHOclass0305)
# 
# 
# dds3 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDEGRADE,
#                                        colData = metadataDEGRADE,
#                                        design= ~Grade)
# dds3 <-DESeq(dds3)
# vsd3 <- DESeq2::vst(dds3)
# readcount.vstDEGRADE <- assay(vsd3) %>% as.data.frame 
# 
# 
# 
# deseq2ResultsGRADE <- results(dds3) %>% 
#   as.data.frame()%>% 
#   rownames_to_column("gene_id") %>% 
#   left_join(gene.annot2)
# 
# rm(readcountsDEGRADE, dds3, vsd3)
# 
# 
# save(deseq2ResultsGRADE, file = "deseq2ResultsGRADE.Rdata")
# 
# 
# deseq2ResultsSIGNGRADE <- deseq2ResultsGRADE%>% 
#   dplyr::filter(padj <0.05) %>% 
#   dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
#   left_join(gene.annot2)
# 
# signupRNAGRADE <- deseq2ResultsSIGNGRADE %>% 
#   dplyr::filter(log2FoldChange >0)
# 
# signdownRNAGRADE <- deseq2ResultsSIGNGRADE %>% 
#   dplyr::filter(log2FoldChange <0)
# 
# 
# 
# 
# 
# upprotRNA3 <- signupRNAGRADE %>% 
#   inner_join(diffprot3signup, by = c("gene_name"= "Protein"))
# 
# downprotRNA3 <- signdownRNAGRADE %>% 
#   inner_join(diffprot3signdown, by = c("gene_name"= "Protein"))
# 
# 
# summary(deseq2Results)
# plotMA(deseq2Results)
# 
# rm(vsd2,tmp, dds2, readcount.vst2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###########PLOTJES
#
DNA_binding <- read_csv("DNA_binding.csv")

library(data.table)

deseq2Results %>% dplyr::filter(contains('gene'))

HOXgenes <- deseq2Results %>% 
  filter(str_detect(gene_name, "HOX")) %>% 
  dplyr::select(gene_name)

histones <- deseq2Results %>% 
  filter(str_detect(gene_name, "^[H][0-9]")) %>% 
  dplyr::select(gene_name)

celldeseq2volcano <- deseq2Results %>%
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>%
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% HOXgenes$gene_name, "HOX genes", colors)) %>%
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & log2FoldChange > 0.75 & gene_name %in% histones$gene_name, "Histones", colors)) %>% 
   dplyr::mutate(size = ifelse(colors == "HOX genes" | colors == "Cell cycle" | colors == "Histones", 1.5,1))



pdf("Volcano_RNA.pdf", width = 10, height = 8)
ggplot(celldeseq2volcano, aes(log2FoldChange, -log10(padj), color = colors, label = gene_name)) +
  geom_point(data = subset(celldeseq2volcano, colors == "Not Significant"), size = 1)+
  geom_point(data = subset(celldeseq2volcano, colors == "Significant"), size = 1)+
  geom_point(data = subset(celldeseq2volcano, colors == "HOX genes"), size = 2)+
  geom_point(data = subset(celldeseq2volcano, colors == "Cell cycle"), size = 2)+
  geom_point(data = subset(celldeseq2volcano, colors == "Histones"), size = 2)+
  scale_color_manual(values = c("gray87", "gray40", "royalblue", "purple", "brown"),
                   name = "colors",
                   breaks=c("Not Significant", "Significant", "HOX genes", "Cell cycle", "Histones"),
                   labels = c("Not Significant", "Significant", "HOX genes", "Cell cycle", "Histones"))+
   theme_bw()+
   geom_text_repel(data=subset(celldeseq2volcano, colors == "Cell cycle" | colors == "HOX genes" | colors == "Histones"), aes(label=gene_name), color = "black", max.overlaps = 5) +
   geom_hline(yintercept=1.3, linetype="dashed")+
   geom_vline(xintercept = 0.75, linetype="dashed")+
   geom_vline(xintercept = -0.75, linetype="dashed")+
   labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between methylation classifier LG (n=112) versus HG (n=41)')+
   xlab('Log fold change')+
   ylab('-log10(adjusted p-value)')+
   xlim(-6.5,6.5)
   

dev.off()


celldeseq2volcano_predictive <- deseq2ResultsLGHG %>%
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>%
   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & 
                                   gene_name %in% c("CENPF","TOP2A", "ASPM","MKI67","FOXM1"), "Predictive cell cycle", colors)) %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & 
                                  gene_name %in% c("MCM10","ADAMTS7", "MNX1", "H3C2", "CRNDE"), "Predictive", colors)) %>% 
   dplyr::mutate(size = ifelse(colors == "Predictive" | colors == "Predictive cell cycle", 3,1))
 
 pdf("Volcano_RNA_predictive.pdf", width = 10, height = 7)
 
 ggplot(celldeseq2volcano_predictive, aes(log2FoldChange, -log10(padj))) +
   geom_point(aes(col=colors, size = size)) +
   scale_radius(range = c(1,2))+
   guides(size = "none")+
   scale_color_manual(values=c("gray87", "red","purple", "gray60")) +
   theme_bw()+
   geom_text_repel(data=subset(celldeseq2volcano_predictive, colors == "Predictive" | colors == "Predictive cell cycle"), aes(label=gene_name), max.overlaps = 5) +
   geom_hline(yintercept=1.3, linetype="dashed")+
   geom_vline(xintercept = 0.75, linetype="dashed")+
   geom_vline(xintercept = -0.75, linetype="dashed")+
   labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between methylation classifier LG (n=77) versus HG (n=34)')+
   xlab('Log fold change')+
   ylab('-log10(adjusted p-value)')+
   xlim(-7,7)
 
 dev.off()
 
 
#
# 
# 
# ###############GRADE
# celldeseq3volcano <- deseq2ResultsGRADE %>% 
#   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
#   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colors)) %>% 
#   dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>% 
#   dplyr::mutate(size = ifelse(colors == "DNA binding" | colors == "Cell cycle", 1.5,1))
# 
# 
# ggplot(celldeseq3volcano, aes(log2FoldChange, -log10(padj))) + 
#   geom_point(aes(col=colors, size = size)) + 
#   scale_radius(range = c(1,2))+ 
#   guides(size = "none")+
#   scale_color_manual(values=c("purple","royalblue", "gray87", "gray40")) + 
#   theme_bw()+
#   geom_text_repel(data=subset(celldeseq3volcano, colors == "Cell cycle" | colors == "DNA binding"), aes(label=gene_name), max.overlaps = 4) +
#   geom_hline(yintercept=1.3, linetype="dashed")+
#   geom_vline(xintercept = 0.75, linetype="dashed")+
#   geom_vline(xintercept = -0.75, linetype="dashed")+
#   labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between WHO2021 grade 2&3 (n=121) vs grade 4 (n=40)')+
#   xlab('Log fold change')+
#   ylab('-log10(adjusted p-value)')+
#   xlim(-7,7)
# 
# 
#enrichment

enrichsignup <- gost(signupRNA$gene_name, organism = "hsapiens", ordered_query = F, significant = T, correction_method = "fdr")

resenrichsignup <- enrichsignup$result


resenrichsignup <- resenrichsignup[order(resenrichsignup$p_value),]

barenrichsignup <- names(resenrichsignup$term_id)
edo <- enrichDGN(barenrichsignup)

gostplot(enrichsignup)


### METHYLATION


upRNA_prenames <- signupRNA %>%
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
upRNA_names <- gconvert(upRNA_prenames$gene_id)

up_gp = gost(list(upRNA_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
up_gp = up_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
up_gp <- up_gp[order(up_gp$p_value),]
up_gp$GeneRatio = paste0(up_gp$intersection_size, "/", up_gp$query_size)
up_gp$BgRatio = paste0(up_gp$term_size, "/", up_gp$effective_domain_size)
names(up_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
up_gp$geneID = gsub(",","/", up_gp$geneID)
row.names(up_gp) = up_gp$ID
up_gp_cluster = new("compareClusterResult", compareClusterResult = up_gp)
up_gp_enrich = new("enrichResult", result = up_gp)
up_gp_enrich@result$p.adjust<- up_gp_enrich@result$p.adjust %>%
  signif(., 2)
result_up_gp <- up_gp_enrich@result
enrichplot::dotplot(up_gp_cluster)
barplot(up_gp_enrich, showCategory = 10, fontsize = 5)


up_gp_gg <-ggbarplot(result_up_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")
pdf("upregRNA.pdf", width = 6, height = 4)
ggpar(up_gp_gg, title = "Upregulated pathways RNA for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count",
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()

rm(up_gp, upRNA_prenames, upRNA_names, up_gp_cluster, up_gp_enrich, up_gp_gg, result_up_gp) 
 
 ##DOWN
downRNA_prenames <- signdownRNA %>% 
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
downRNA_names <- gconvert((downRNA_prenames$gene_id))
 
down_gp = gost(list(downRNA_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
down_gp = down_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
down_gp <- down_gp[order(down_gp$p_value),]
down_gp$GeneRatio = paste0(down_gp$intersection_size, "/", down_gp$query_size)
down_gp$BgRatio = paste0(down_gp$term_size, "/", down_gp$effective_domain_size)
names(down_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
down_gp$geneID = gsub(",","/", down_gp$geneID)
row.names(down_gp) = down_gp$ID
down_gp_cluster = new("compareClusterResult", compareClusterResult = down_gp)
down_gp_enrich = new("enrichResult", result = down_gp)
down_gp_enrich@result$p.adjust<- down_gp_enrich@result$p.adjust %>% 
 signif(., 2)
result_down_gp <- down_gp_enrich@result

down_gp_gg <-ggbarplot(result_down_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")

 
pdf("downregRNA.pdf", width = 6, height = 4)

ggpar(down_gp_gg, title = "Downregulated pathways RNA for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
       font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)
 

dev.off()
 
#ggplot(result_down_gp[1:10,], aes (x = Count, y = Description, fill = p.adjust))+
#  geom_col()+
#  scale_color_gradient2(low = "black", high= "white")
# 
# 
# save(result_down_gp, file = "result_down_gp.Rdata")
# save(result_up_gp, file = "result_up_gp.Rdata")
# 
# 
# 

# metadataDE2 <- metadataDE %>% 
#   dplyr::filter(Sample_Name %in% c("210_R2", "177_R2", "023_R2") ==F)
# 
# 
# readcountsDE2 <- readcountsDE %>% 
#   dplyr::mutate("104059-003-035" = NULL) %>% 
#   dplyr::mutate("104059-002-187" = NULL) %>% 
#   dplyr::mutate("104059-002-037" = NULL)
# 
# 
# dds5 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE2,
#                                        colData = metadataDE2,
#                                        design= ~methylation.sub.diagnosis + Sample_Type)
# dds5 <-DESeq(dds5)
# 
# coefdds5 <- coef(dds5)
# 
# coefssignupprimrecur <- coefdds5 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("gene_id") %>% 
#   inner_join(deseq2ResultsSIGNpr)
# 
# coefssignupprimrecur <- coefssignupprimrecur[order(coefssignupprimrecur$log2FoldChange),]
# 
# coefssignupprimrecur <- coefssignupprimrecur[match(o, coefssignupprimrecur$gene_id),]
# 
# 
# barplot((coefssignupprimrecur$methylation.sub.diagnosis_A_IDH_HG_vs_A_IDH/3), ylim = c(-4,4), ylab = "High vs Low grade", cex.lab =1.3)
# 
# barplot((coefssignupprimrecur$Sample_Type_R_vs_I/3), ylim = c(-4,4), ylab = "Primary vs Recurrent", cex.lab =1.3)
# 
# 
# 
# 
# 
# corplotprimrecur <- readcount.vst %>% 
#   rownames_to_column("gene_id")
# 
# primrecurRNA <- coefssignupprimrecur %>% 
#   dplyr::select(gene_id) %>% 
#   left_join(corplotprimrecur) %>% 
#   tibble::column_to_rownames("gene_id") %>% 
#   as.data.frame()
# 
# 
# 
# recursiveCorPlot(primrecurRNA, metadata, font_scale = 2, legend_scale =3, method="ward.D2", return_h_object = FALSE, caption=NULL)
# 
# 
# plt <- primrecurRNA %>%
#   as.matrix %>%
#   t() %>%
#   stats::cor()
# 
# h <- stats::hclust(stats::as.dist(1 - stats::cor(plt)), method="ward.D2" ) # recursive cor-based cluastering !!!
# #h <- stats::hclust( stats::as.dist(1 - plt) , method = method ) # regular cor-based clustering
# 
# o <- h$labels[h$order] %>% rev()
# 
# 
# 
# 
# 
# dds6 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE,
#                                        colData = metadataDE,
#                                        design= ~Grade + Sample_Type)
# dds6 <-DESeq(dds6)
# 
# coefdds6 <- coef(dds6)
# 
# 
# 
# 
# readcount.vstDE <- assay(vsd2) %>% as.data.frame 
# 
# 
# signupRNAboth <- inner_join(signupRNAGRADE, signupRNA, by = "gene_name")
# 
# 
# 
# #----------------------------------------------------


metadata_DE_paired <- metadataDE %>%
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Sample_Name)) %>%
  dplyr::mutate(batch = ifelse(Sample_Number %in% c("207", "054", "048","173","168","158","143","141","136","135","134","130","128","126","121","118","117","116",
                                                            "113","111","108","107","104","049","039","035","033","030","025","020","015","010","002"),"001", Sample_Number))%>%
  dplyr::mutate(batch = paste0("s",batch)) %>%
  dplyr::mutate(batch = factor(batch)) %>% 
  dplyr::mutate(methylation.sub.diagnosis  = factor(methylation.sub.diagnosis , levels=c('A_IDH','A_IDH_HG')))


stopifnot(length(metadata_DE_paired$Sample_Number[duplicated(metadata_DE_paired$Sample_Number)]) -sum(metadata_DE_paired$Sample_Number[duplicated(metadata_DE_paired$Sample_Number)] %>% duplicated()) +1 == length(metadata_DE_paired$batch %>% levels))

dds_paired <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE,
                                             colData = metadata_DE_paired,
                                             design= ~batch + methylation.sub.diagnosis)

dds_paired <-DESeq(dds_paired)


deseq2Results_paired <- results(dds_paired) %>%
  as.data.frame()%>%
  rownames_to_column("gene_id") %>%
  left_join(gene.annot2)

# 
# 
# 
# ##DESEQ OP PATIENTEN DIE IN PROTEOME VOORKOMEN
# readcountspatprot <- readcountsDE %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("GS_ID") %>% 
#   dplyr::filter(GS_ID %in% sortpatientid) %>% 
#   column_to_rownames("GS_ID") %>% 
#   t() %>% 
#   as.data.frame()
# 
# metadatapatprot <- metadataDE %>% 
#   dplyr::filter(GS_ID %in% sortpatientid)
# 
# 
# stopifnot(length(colnames(readcountspatprot)) == length(metadatapatprot$GS_ID))
# 
# ddspatprot <- DESeq2::DESeqDataSetFromMatrix(countData = readcountspatprot,
#                                        colData = metadatapatprot,
#                                        design= ~methylation.sub.diagnosis)
# ddspatprot <-DESeq(ddspatprot)
# vsdpatprot <- DESeq2::vst(ddspatprot)
# readcount.vstpatprot <- assay(vsdpatprot) %>% as.data.frame
# 
# deseq2Resultspatprot <- results(ddspatprot) %>%
#   as.data.frame()%>%
#   rownames_to_column("gene_id") %>%
#   left_join(gene.annot2)

deseqpaired <- res_deseq %>% 
  as.data.frame()%>%
  rownames_to_column("gene_id") %>%
  left_join(gene.annot2)

deseq2ResultsSIGN_paired <- deseqpaired%>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  left_join(gene.annot2)

signupRNA_paired <- deseq2ResultsSIGN_paired %>% 
  dplyr::filter(log2FoldChange >0)

signdownRNA_paired <- deseq2ResultsSIGN_paired %>% 
  dplyr::filter(log2FoldChange <0)


##ONLY SAME PATIENTS IN LG VS HG
metadataLGHG<- metadataDE %>% 
  dplyr::select(GLASS_ID, methylation.sub.diagnosis, GS_ID) %>% 
  dplyr::filter(!(GS_ID == "104059-002-007")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-008")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-009")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-014")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-019")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-021")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-030")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-038")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-035")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-044")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-049")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-052")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-055")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-060")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-067")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-072")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-077")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-104")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-108")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-111")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-124")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-127")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-135")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-137")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-152")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-156")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-158")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-162")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-171")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-174")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-178")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-181")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-187")) %>% 
  dplyr::filter(!(GS_ID == "104059-002-188")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-004")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-006")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-012")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-026")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-030")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-035")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-040")) %>% 
  dplyr::filter(!(GS_ID == "104059-003-046")) %>% 
  dplyr::mutate(Sample_Number = str_sub(GLASS_ID, start = -3)) %>% 
  dplyr::mutate(Sample_Number = ifelse(Sample_Number %!in% c("017", "023", "029", "032", "034", "043", "045", "125", "129", "142", "145", "147",
                                                             "149", "156", "163", "171", "172", "174", "044", "115", "146", "202", "211"), Sample_Number, "001"))


  

metadataLGHG$GLASS_ID[duplicated(metadataLGHG$GLASS_ID)]

readcountsLGHG <- readcountsDE %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% metadataLGHG$GS_ID) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()


ddsLGHG <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsLGHG,
                                      colData = metadataLGHG,
                                      design= ~methylation.sub.diagnosis)
ddsLGHG <-DESeq(ddsLGHG)
vsdLGHG <- DESeq2::vst(ddsLGHG)
readcount.vstLGHG <- assay(vsdLGHG) %>% as.data.frame 


deseq2ResultsLGHG <- results(ddsLGHG) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)


signupLGHG <- deseq2ResultsLGHG %>% 
  dplyr::filter(padj <0.05 & log2FoldChange >0.75)

signdownLGHG <- deseq2ResultsLGHG %>% 
  dplyr::filter(padj <0.05 & log2FoldChange < -0.75)

celldeseq2volcanoLGHG <- deseq2ResultsLGHG %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% HOXgenes$gene_name, "HOX genes", colors)) %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>%
  dplyr::mutate(colors =ifelse(padj <0.05 & log2FoldChange > 0.75 & gene_name %in% histones$gene_name, "Histones", colors)) %>% 
  dplyr::mutate(colors =ifelse(padj <0.05 & log2FoldChange < -0.75 & gene_name %in% ion_transport$name, "Ion transport", colors)) %>% 
  dplyr::mutate(size = ifelse(colors == "HOX genes" | colors == "Cell cycle" | colors == "Histones" | colors == "Ion transport", 1.5,1))



pdf("Volcano_RNA.pdf", width = 10, height = 8)
ggplot(celldeseq2volcanoLGHG, aes(log2FoldChange, -log10(padj), color = colors, label = gene_name)) +
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "Not Significant"), size = 1)+
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "Significant"), size = 1)+
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "HOX genes"), size = 2)+
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "Cell cycle"), size = 2)+
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "Histones"), size = 2)+
  geom_point(data = subset(celldeseq2volcanoLGHG, colors == "Ion transport"), size = 2)+
  scale_color_manual(values = c("gray87", "gray40", "royalblue", "purple", "brown", "red"),
                     name = "colors",
                     breaks=c("Not Significant", "Significant", "HOX genes", "Cell cycle", "Histones", "Ion transport"),
                     labels = c("Not Significant", "Significant", "HOX genes", "Cell cycle", "Histones", "Ion transport"))+
  theme_bw()+
  geom_text_repel(data=subset(celldeseq2volcanoLGHG, colors == "Cell cycle" | colors == "HOX genes" | colors == "Histones"| colors == "Ion transport"), aes(label=gene_name), color = "black", max.overlaps = 5) +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between methylation classifier LG (n=77) versus HG (n=34)')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-6.5,6.5)
dev.off()

##ENRICHMENT

upRNA_prenames <- signupLGHG %>%
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
upRNA_names <- gconvert(upRNA_prenames$gene_id)

up_gp = gost(list(upRNA_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
up_gp = up_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
up_gp <- up_gp[order(up_gp$p_value),]
up_gp$GeneRatio = paste0(up_gp$intersection_size, "/", up_gp$query_size)
up_gp$BgRatio = paste0(up_gp$term_size, "/", up_gp$effective_domain_size)
names(up_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
up_gp$geneID = gsub(",","/", up_gp$geneID)
row.names(up_gp) = up_gp$ID
up_gp_cluster = new("compareClusterResult", compareClusterResult = up_gp)
up_gp_enrich = new("enrichResult", result = up_gp)
up_gp_enrich@result$p.adjust<- up_gp_enrich@result$p.adjust %>%
  signif(., 2)
result_up_gp <- up_gp_enrich@result
enrichplot::dotplot(up_gp_cluster)
barplot(up_gp_enrich, showCategory = 10, fontsize = 5)


up_gp_gg <-ggbarplot(result_up_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")
pdf("upregRNA.pdf", width = 6, height = 4)
ggpar(up_gp_gg, title = "Upregulated pathways RNA for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count",
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()

rm(up_gp, upRNA_prenames, upRNA_names, up_gp_cluster, up_gp_enrich, up_gp_gg, result_up_gp) 

##DOWN
downRNA_prenames <- signdownLGHG %>% 
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
downRNA_names <- gconvert((downRNA_prenames$gene_id))

down_gp = gost(list(downRNA_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
down_gp = down_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
down_gp <- down_gp[order(down_gp$p_value),]
down_gp$GeneRatio = paste0(down_gp$intersection_size, "/", down_gp$query_size)
down_gp$BgRatio = paste0(down_gp$term_size, "/", down_gp$effective_domain_size)
names(down_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
down_gp$geneID = gsub(",","/", down_gp$geneID)
row.names(down_gp) = down_gp$ID
down_gp_cluster = new("compareClusterResult", compareClusterResult = down_gp)
down_gp_enrich = new("enrichResult", result = down_gp)
down_gp_enrich@result$p.adjust<- down_gp_enrich@result$p.adjust %>% 
  signif(., 2)
result_down_gp <- down_gp_enrich@result

down_gp_gg <-ggbarplot(result_down_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")


pdf("downregRNA.pdf", width = 6, height = 4)

ggpar(down_gp_gg, title = "Downregulated pathways RNA for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)


dev.off()



signupboth <- signupLGHG %>% 
  inner_join(diffprot2signupLGHG, by = c("gene_name" = "Protein"))

signdownboth <- signdownLGHG %>% 
  inner_join(diffprot2signdownLGHG, by = c("gene_name" = "Protein"))
  

upprimrecurLGHG <- signupRNAprimrecur_lastres %>% 
  inner_join(signupLGHG, by = "gene_id")

downprimrecurLGHG <- signdownRNAprimrecur_lastres %>% 
  inner_join(signdownLGHG, by = "gene_id")


##

ddsLGHG_paired <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsLGHG,
                                          colData = metadataLGHG,
                                          design= ~Sample_Number + methylation.sub.diagnosis)
ddsLGHG_paired <-DESeq(ddsLGHG_paired)
vsdLGHG_paired <- DESeq2::vst(ddsLGHG_paired)
readcount.vstLGHG_paired <- assay(vsdLGHG_paired) %>% as.data.frame 


RNA  <- c(1699, 1052)
prot <- c(82, 60)
both <- c(15,8)

chidf <- data.frame(RNA, prot, both, row.names = c("up", "down"))

chisq <- chisq.test(chidf)
chisq$expected
chisq$p.value
