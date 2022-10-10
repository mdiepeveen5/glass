library(tidyverse)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

#---- load data
library(readODS)
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
  dplyr::filter(GS_ID %in% metadataDE$GS_ID) %>% 
  tibble::column_to_rownames('GS_ID') %>% 
  t() %>% 
  as.data.frame()


#  dplyr::filter(GS_ID %in% metadataDE$GS_ID) %>% 

load(file = "metadataALL.Rdata")


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
  dplyr::filter(methylation.sub.diagnosis %in% c("A_IDH_HG", "A_IDH")) %>% 
  left_join(WHOclass0305)


#maybe leave inner_join out if doesn't work 

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
#firsts run readcountsDE without GS_ID %in% metadataDE, then run metadataDE, then run readcountsDE with GS_ID %in% metadataDE

WHOclass0305<- read.csv("WHOclassification_03052022.csv") %>% 
  dplyr::filter(Sample_Name %in% c("007_R4","026_R2", "219_R3", "216_R3") ==F) %>% 
  dplyr::mutate(Grade = ifelse(WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 4", "Grade4", "Grade2or3"))



###GRADE
readcountsDEGRADE <- tmp %>% 
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
  dplyr::filter(GS_ID %in% metadataDEGRADE$GS_ID) %>% 
  tibble::column_to_rownames('GS_ID') %>% 
  t() %>% 
  as.data.frame()

metadataDEGRADE<-
  data.frame(GS_ID = colnames(readcountsDEGRADE)) %>% 
  dplyr::left_join(
    read.csv("~/projects/glass/data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") ,
    by=c('GS_ID'='GS_ID')
  ) %>% 
  dplyr::mutate(Exclude=NULL) %>%   
  dplyr::mutate(Surgery_ID=NULL) %>% 
  dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex)) %>% 
  dplyr::left_join(metadataALL) %>% 
  dplyr::inner_join(WHOclass0305)


dds3 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDEGRADE,
                                       colData = metadataDEGRADE,
                                       design= ~Grade)
dds3 <-DESeq(dds3)
vsd3 <- DESeq2::vst(dds3)
readcount.vstDEGRADE <- assay(vsd3) %>% as.data.frame 



deseq2ResultsGRADE <- results(dds3) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)

rm(readcountsDEGRADE, dds3, vsd3)


save(deseq2ResultsGRADE, file = "deseq2ResultsGRADE.Rdata")


deseq2ResultsSIGNGRADE <- deseq2ResultsGRADE%>% 
  dplyr::filter(padj <0.05) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.75) %>% 
  left_join(gene.annot2)

signupRNAGRADE <- deseq2ResultsSIGNGRADE %>% 
  dplyr::filter(log2FoldChange >0)

signdownRNAGRADE <- deseq2ResultsSIGNGRADE %>% 
  dplyr::filter(log2FoldChange <0)





upprotRNA3 <- signupRNAGRADE %>% 
  inner_join(diffprot3signup, by = c("gene_name"= "Protein"))

downprotRNA3 <- signdownRNAGRADE %>% 
  inner_join(diffprot3signdown, by = c("gene_name"= "Protein"))


summary(deseq2Results)
plotMA(deseq2Results)

rm(vsd2,tmp, dds2, readcount.vst2)









###########PLOTJES
cell_cycle <- read_csv("cell_cycle.csv")

DNA_binding <- read_csv("DNA_binding.csv")

cellcycledeseq2 <- inner_join(as.data.frame(deseq2Results$gene_name), as.data.frame(cell_cycle$name), by = c("deseq2Results$gene_name" = "cell_cycle$name")) %>% 
  dplyr::rename("gene" = "deseq2Results$gene_name")

DNAbindingdeseq2 <- inner_join(as.data.frame(deseq2Results$gene_name), as.data.frame(DNA_binding$name), by = c("deseq2Results$gene_name" = "DNA_binding$name")) %>% 
  dplyr::rename("gene" = "deseq2Results$gene_name")

celldeseq2volcano <- deseq2Results

celldeseq2volcano$group[celldeseq2volcano$gene_name %in% cycling.cell.markers$gene_name] <- "Cell cycle"
celldeseq2volcano$group[celldeseq2volcano$gene_name %in% DNAbindingdeseq2$gene] <- "DNA binding"


keyvals <- ifelse(
  celldeseq2volcano$group == "DNA binding" & deseq2Results$log2FoldChange >abs(0.75), 'royalblue',
  ifelse(celldeseq2volcano$group == "Cell cycle" & deseq2Results$log2FoldChange > abs(0.75), "purple", "grey"))
keyvals[is.na(keyvals)] <- "grey"
names(keyvals)[keyvals == 'royalblue'] <- 'DNA binding'
names(keyvals)[keyvals == 'purple'] <- 'Cell cycle'
names(keyvals)[keyvals == "grey"] <- 'Other'

deseq2Results$gene_name[which(names(keyvals) %in% c("DNA binding", "Cell cycle"))]


logfcup2 <- celldeseq2volcano %>% 
  dplyr::filter(!is.na(group)) %>% 
  dplyr::filter(log2FoldChange >2) %>% 
  dplyr::select(gene_name)


#ADJUSTED P VAL 
EnhancedVolcano(celldeseq2volcano,
                lab = celldeseq2volcano$gene_name,
                selectLab = logfcup2$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Differentially expressed genes for RNA-sequencing data',
                colCustom = keyvals,
                pointSize = c(ifelse(celldeseq2volcano$log2FoldChange>2, 3, 1)),
                pCutoff = 0.05,
                FCcutoff = 0.75,
                labSize = 6.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colAlpha = 1,
                drawConnectors = T)
library(EnhancedVolcano)

celldeseq2volcano <- deseq2Results %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colours)) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colours)) %>% 
  dplyr::mutate(size = ifelse(colours == "DNA binding" | colours == "Cell cycle", 1.5,1))


ggplot(celldeseq2volcano, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col=colours, size = size)) + 
  scale_radius(range = c(1,2))+ 
  guides(size = "none")+
  scale_color_manual(values=c("purple","royalblue", "gray87", "gray40")) + 
  theme_bw()+
  geom_text_repel(data=subset(celldeseq2volcano, colours == "Cell cycle" | colours == "DNA binding"), aes(label=gene_name), max.overlaps = 4) +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between methylation classifier LG (n=112) vs HG (n=41)')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-7,7)



###############GRADE
celldeseq3volcano <- deseq2ResultsGRADE %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75,  "Significant", "Not Significant")) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% DNA_binding$name, "DNA binding", colours)) %>% 
  dplyr::mutate(colours =ifelse(padj <0.05 & abs(log2FoldChange) > 0.75 & gene_name %in% cycling.cell.markers$gene_name, "Cell cycle", colours)) %>% 
  dplyr::mutate(size = ifelse(colours == "DNA binding" | colours == "Cell cycle", 1.5,1))


ggplot(celldeseq3volcano, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col=colours, size = size)) + 
  scale_radius(range = c(1,2))+ 
  guides(size = "none")+
  scale_color_manual(values=c("purple","royalblue", "gray87", "gray40")) + 
  theme_bw()+
  geom_text_repel(data=subset(celldeseq3volcano, colours == "Cell cycle" | colours == "DNA binding"), aes(label=gene_name), max.overlaps = 4) +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for RNA data', subtitle = 'Between WHO2021 grade 2&3 (n=121) vs grade 4 (n=40)')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-7,7)


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

ggpar(up_gp_gg, title = "Upregulated pathways RNA between methylation classifier LG vs HG", font.main = c(23, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(20, "plain", "black"), font.y = c(20, "plain", "black"), font.ytickslab = c(16, "plain", "black"))



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


ggpar(down_gp_gg, title = "Downregulated pathways RNA between methylation classifier LG vs HG", font.main = c(23, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(20, "plain", "black"), font.y = c(20, "plain", "black"), font.ytickslab = c(16, "plain", "black"))



ggplot(result_down_gp[1:10,], aes (x = Count, y = Description, fill = p.adjust))+
  geom_col()+
  scale_color_gradient2(low = "black", high= "white")


save(result_down_gp, file = "result_down_gp.Rdata")
save(result_up_gp, file = "result_up_gp.Rdata")



####################GRADE

upRNAGRADE_prenames <- signupRNAGRADE %>% 
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
upRNA_namesGRADE <- gconvert(upRNAGRADE_prenames$gene_id)

up_gpGRADE = gost(list(upRNA_namesGRADE$name), multi_query = F, evcodes = T, sources = "GO:BP")
up_gpGRADE = up_gpGRADE$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
up_gpGRADE <- up_gpGRADE[order(up_gpGRADE$p_value),]
up_gpGRADE$GeneRatio = paste0(up_gpGRADE$intersection_size, "/", up_gpGRADE$query_size)
up_gpGRADE$BgRatio = paste0(up_gpGRADE$term_size, "/", up_gpGRADE$effective_domain_size)
names(up_gpGRADE) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
up_gpGRADE$geneID = gsub(",","/", up_gpGRADE$geneID)
row.names(up_gpGRADE) = up_gpGRADE$ID
up_gpGRADE_cluster = new("compareClusterResult", compareClusterResult = up_gpGRADE)
up_gpGRADE_enrich = new("enrichResult", result = up_gpGRADE)
up_gpGRADE_enrich@result$p.adjust<- up_gpGRADE_enrich@result$p.adjust %>% 
  signif(., 2)
result_up_gpGRADE <- up_gpGRADE_enrich@result
enrichplot::dotplot(up_gpGRADE_cluster)
barplot(up_gpGRADE_enrich, showCategory = 10, fontsize = 5)

up_gpGRADE_gg <-ggbarplot(result_up_gpGRADE[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")

ggpar(up_gpGRADE_gg, title = "Upregulated pathways RNA between WHO2021 grade 2&3 vs grade 4", font.main = c(18, "bold", "#0d6d42"), xlab = "Description", ylab= "Count", 
      font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"))



downRNAGRADE_prenames <- signdownRNAGRADE %>% 
  dplyr::mutate(gene_id = gsub ("\\..*", "", gene_id))
downRNA_namesGRADE <- gconvert((downRNAGRADE_prenames$gene_id))

down_gpGRADE = gost(list(downRNA_namesGRADE$name), multi_query = F, evcodes = T, sources = "GO:BP")
down_gpGRADE = down_gpGRADE$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
down_gpGRADE <- down_gpGRADE[order(down_gpGRADE$p_value),]
down_gpGRADE$GeneRatio = paste0(down_gpGRADE$intersection_size, "/", down_gpGRADE$query_size)
down_gpGRADE$BgRatio = paste0(down_gpGRADE$term_size, "/", down_gpGRADE$effective_domain_size)
names(down_gpGRADE) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
down_gpGRADE$geneID = gsub(",","/", down_gpGRADE$geneID)
row.names(down_gpGRADE) = down_gpGRADE$ID
down_gpGRADE_cluster = new("compareClusterResult", compareClusterResult = down_gpGRADE)
down_gpGRADE_enrich = new("enrichResult", result = down_gpGRADE)
down_gpGRADE_enrich@result$p.adjust<- down_gpGRADE_enrich@result$p.adjust %>% 
  signif(., 2)
result_down_gpGRADE <- down_gpGRADE_enrich@result

down_gpGRADE_gg <-ggbarplot(result_down_gpGRADE[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")


ggpar(down_gpGRADE_gg, title = "Downregulated pathways RNA between WHO2021 grade 2&3 vs grade 4", font.main = c(18, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"))



metadataDE2 <- metadataDE %>% 
  dplyr::filter(Sample_Name %in% c("210_R2", "177_R2", "023_R2") ==F)


readcountsDE2 <- readcountsDE %>% 
  dplyr::mutate("104059-003-035" = NULL) %>% 
  dplyr::mutate("104059-002-187" = NULL) %>% 
  dplyr::mutate("104059-002-037" = NULL)


dds5 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE2,
                                       colData = metadataDE2,
                                       design= ~methylation.sub.diagnosis + Sample_Type)
dds5 <-DESeq(dds5)

coefdds5 <- coef(dds5)

coefssignupprimrecur <- coefdds5 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  inner_join(deseq2ResultsSIGNpr)

coefssignupprimrecur <- coefssignupprimrecur[order(coefssignupprimrecur$log2FoldChange),]

coefssignupprimrecur <- coefssignupprimrecur[match(o, coefssignupprimrecur$gene_id),]


barplot((coefssignupprimrecur$methylation.sub.diagnosis_A_IDH_HG_vs_A_IDH/3), ylim = c(-4,4), ylab = "High vs Low grade", cex.lab =1.3)

barplot((coefssignupprimrecur$Sample_Type_R_vs_I/3), ylim = c(-4,4), ylab = "Primary vs Recurrent", cex.lab =1.3)





corplotprimrecur <- readcount.vst %>% 
  rownames_to_column("gene_id")

primrecurRNA <- coefssignupprimrecur %>% 
  dplyr::select(gene_id) %>% 
  left_join(corplotprimrecur) %>% 
  tibble::column_to_rownames("gene_id") %>% 
  as.data.frame()



recursiveCorPlot(primrecurRNA, metadata, font_scale = 2, legend_scale =3, method="ward.D2", return_h_object = FALSE, caption=NULL)


plt <- primrecurRNA %>%
  as.matrix %>%
  t() %>%
  stats::cor()

h <- stats::hclust(stats::as.dist(1 - stats::cor(plt)), method="ward.D2" ) # recursive cor-based cluastering !!!
#h <- stats::hclust( stats::as.dist(1 - plt) , method = method ) # regular cor-based clustering

o <- h$labels[h$order] %>% rev()





dds6 <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE,
                                       colData = metadataDE,
                                       design= ~Grade + Sample_Type)
dds6 <-DESeq(dds6)

coefdds6 <- coef(dds6)




readcount.vstDE <- assay(vsd2) %>% as.data.frame 


signupRNAboth <- inner_join(signupRNAGRADE, signupRNA, by = "gene_name")



#----------------------------------------------------


metadata_DE_paired <- metadataDE %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Sample_Name)) %>% 
  dplyr::mutate(Sample_Number = ifelse(Sample_Number %in% c("207", "054", "048","173","168","158","143","141","136","135","134","130","128","126","121","118","117","116",
                                                            "113","111","108","107","104","049","039","035","033","030","025","020","015","010","002"),"001", Sample_Number))%>%
  dplyr::mutate(batch = paste0("s",Sample_Number)) %>%
  dplyr::mutate(methylation.sub.diagnosis  = factor(methylation.sub.diagnosis , levels=c('A_IDH','A_IDH_HG')))



dds_paired <- DESeq2::DESeqDataSetFromMatrix(countData = readcountsDE,
                                             colData = metadata_DE_paired,
                                             design= ~batch + methylation.sub.diagnosis)

dds_paired <-DESeq(dds_paired,parallel=F)


deseq2Results_paired <- results(dds_paired) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)




##DESEQ OP PATIENTEN DIE IN PROTEOME VOORKOMEN
readcountspatprot <- readcountsDE %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% sortpatientid) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()

metadatapatprot <- metadataDE %>% 
  dplyr::filter(GS_ID %in% sortpatientid)


stopifnot(length(colnames(readcountspatprot)) == length(metadatapatprot$GS_ID))

ddspatprot <- DESeq2::DESeqDataSetFromMatrix(countData = readcountspatprot,
                                       colData = metadatapatprot,
                                       design= ~methylation.sub.diagnosis)
ddspatprot <-DESeq(ddspatprot)
vsdpatprot <- DESeq2::vst(ddspatprot)
readcount.vstpatprot <- assay(vsdpatprot) %>% as.data.frame 

deseq2Resultspatprot <- results(ddspatprot) %>% 
  as.data.frame()%>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2)
