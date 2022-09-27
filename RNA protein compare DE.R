library("ggplot2")
library("ggrepel")
#needs diffprot2 and deseq2Results and gene.annot2

load(file = "diffprot2.Rdata")
load(file = "deseq2Results.Rdata")

diffprotdeseq <- deseq2Results %>% 
  left_join(gene.annot2) %>% 
  inner_join(diffprot2, by = c("gene_name"="Protein"))


diffprotdeseq2 <- deseq2Results %>% 
  left_join(gene.annot2) %>% 
  inner_join(diffprot2, by = c("gene_name"="Protein"))

library(readr)
cell_cycle <- read_csv("cell_cycle.csv")



cellcycleindiffprot <- inner_join(as.data.frame(diffprotdeseq2$gene_name), as.data.frame(cell_cycle$name), by = c("diffprotdeseq2$gene_name" = "cell_cycle$name")) %>% 
  dplyr::rename("gene" = "diffprotdeseq2$gene_name")

save(cellcycleindiffprot, file = "cellcycleindiffprot.Rdata")

diffprotdeseq2$group[diffprotdeseq2$gene_name %in% cell_cycle$name] <- "Cell cycle"

print(as.data.frame(cellcycleindiffprot), row.names = F)
  
diffprotdeseq2$group<- "Not differentially expressed"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c ("RCC2","SRRM1","HMGN2","RCC1","RBBP4","ANP32E","EEF1AKNMT","CDC73","H3-3A","PRPF40A","CNOT9","COL6A3","ACTL6A","SOX2","RFC4","PDS5A","FIP1L1","PPAT","SMARCA5","TGFBI",   
                                                    "TXNDC5","DEK","H1-5","MDC1","MAPK14","MCM3","RPF2","NUP43","RPA3","CBX3","ASNS","BUD31","MCM7","ENY2","RAD21","NCBP1","PRPF4","GOLGA2","GTF3C5","DDX50",   
                                                    "SSRP1","MTA2","SERPINH1","MRE11","DDX47","SMARCD1","SNRPF","TMPO","NOVA1","RSL1D1","MYBBP1A","SUPT6H","ZNF207","LIG3","SMARCE1","VPS25","PSME3","KPNA2","JPT1","HDGFL2",  
                                                    "RAVER1","SMARCA4","NFIX","SAMD1","UBA2","BAX","PCNA","XRN2","RBM12","RPRD1B","COL18A1","COL6A2","PRDX4","USP11","GPKOW","TCEAL4","HMGB3","BGN","HCFC1" )] <- "Differential up"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% cell_cycle$name] <- "Cell cycle"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("ATP1A1","PBXIP1","TNR","SLC1A4","SFXN5","TUBA4A","PLCD1","LSAMP","RAB6B","PCCB","MRAS","NCEH1","PEX5L","MAP6D1","ATP5ME","SLC4A4","ENPP6","SLC25A4","TPPP","CPLANE1",  
                                                   "PURA","KCTD16","CUTA","CPNE5","RALA","GNAI1","ATP5MF","AGBL3","KBTBD11","NEFM","NEFL","ATP6V1H","SYBU","GLIPR2","NIPSNAP3A","PTGDS","CAMK2G","GLUD1","INA","CNTN1",    
                                                   "ATP5F1B","ALDH2","PEBP1","ALDH6A1","NRXN3","TEDC1","CDIP1","ABAT","GNAO1","PLLP","WDR59","CDH13","ALDOC","CNP","CNTNAP1","MBP","TUBB4A","MAG","ATP1A3","SNTA1",    
                                                   "GNAZ","NEFH","PLP1","RAP2C","MT-ATP6"  )] <- "Differential down"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("PCNA", "MCM3","MCM7", "RCC1", "MYBBP1A" )] <- "Paralogs upgerulated in RNA"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("H4C1", "H3-3A", "H2BC17", "H2BC12", "H2AZ1","H2AC4", "H2AC21", "H2AC20", "H2AC11", "H1-5", "H1-4", "H1-3","H1-2", "H1-10","H1-0")] <- "Histones"



ggplot(diffprotdeseq2, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq2, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq2, group == "Differential up"), size=1)+
  geom_point(data = subset(diffprotdeseq2, group == "Differential down"), size=1)+
  geom_point(data = subset(diffprotdeseq2, group == "Histones"), size=3)+
  geom_point(data = subset(diffprotdeseq2, group == "Cell cycle"), size=3)+
  geom_point(data = subset(diffprotdeseq2, group == "Paralogs upgerulated in RNA"), size=3)+
  geom_text_repel(data = subset(diffprotdeseq2, group == "Paralogs upgerulated in RNA"), color = "black")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_abline(intercept = 0, slope = 0.5)+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")+
  scale_color_manual(values = c("grey", "blue", "red", "brown", "mediumpurple", "yellow"),
                     name = "Group",
                     breaks=c("Not differentially expressed", "Differential up", "Differential down", "Histones", "Cell cycle",
                              "Paralogs upgerulated in RNA"),
                     labels = c("Not differentially expressed", "Differential up", "Differential down", "Histones", "Cell cycle",
                                "Paralogs upgerulated in RNA"))

##################
















diffprotdeseq$group<- "Not differentially expressed"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c("VWA1","HSPG2","RCC1","CTPS1","PTGFRN","S100A6","TAGLN2","LAMC1","CHI3L1","NID1","TP53I3","IGKC","PDK1","FN1","COL6A3","GAP43","RFC4","LMNB1","TGFBI","PDLIM7","F13A1","H4C1",
                                                   "H2AC4","H1-3","H2AC11","H2BC12","H1-5","H2BC17","TUBB","CLIC1","HMGA1","MCM3","FABP7","ABRACL","AKAP12","SOD2","FAM126A","HSPB1","COL1A2","VGF","NAMPT","CALU","CALD1","PDIA4",
                                                   "FABP5","FABP4","ANXA1","VIM","ITGB1","CDK1","LDHA","SERPING1" ,"PGM2L1","SERPINH1" ,"CAPN5","MCAM","SLC2A3","A2M","TSFM","LYZ", "LUM","TMPO","COL4A2","ACTN1","SERPINA3","AHNAK2",
                                                    "TEDC1","IGHG3","IGHM","THBS1","ANXA2","ACAN","THOC6","CALB2","SERPINF1","LRRC75A" , "CAVIN1","COL1A1","KPNA2","JPT1","TUBB6","LMNB2","DNAJB1","TPM4","LSM4","FCGBP","SNRPB","PCNA",
                                                   "MMP9","RCAN1","CBR3","COL18A1","COL6A2","LGALS1","FBLN1","PRDX4","BGN","L1CAM" )] <- "UP - RNA"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c("GSTM2","GSTM3","AHCYL1","ATP1A2","TNR","AGT","SLC1A4","CHDH","FAM107A","CADM2","OSBPL11","ALDH1L1","CD38","SLC4A4","ANXA3","CPE","ACSL6","DAAM2","OGFRL1","GJA1","MPP6","ADCYAP1R1",
                                                   "EGFR","EPHX2","ADHFE1","RIDA","TJP2","ALDH1A1","STOM","NEBL","SLC1A2","HEPACAM","METTL7A","ALDH2","AMER2","WASF3","CAB39L","NDRG2","RGS6","ALDH6A1","ALDOC","PPP1R1B","ITGB4","SNTA1",
                                                   "ATP9A","SLC25A18","BCR","MAOA")] <- "DOWN - RNA"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c ("RCC2","SRRM1","HMGN2","RCC1","RBBP4","ANP32E","EEF1AKNMT","CDC73","H3-3A","PRPF40A","CNOT9","COL6A3","ACTL6A","SOX2","RFC4","PDS5A","FIP1L1","PPAT","SMARCA5","TGFBI",   
                                                    "TXNDC5","DEK","H1-5","MDC1","MAPK14","MCM3","RPF2","NUP43","RPA3","CBX3","ASNS","BUD31","MCM7","ENY2","RAD21","NCBP1","PRPF4","GOLGA2","GTF3C5","DDX50",   
                                                    "SSRP1","MTA2","SERPINH1","MRE11","DDX47","SMARCD1","SNRPF","TMPO","NOVA1","RSL1D1","MYBBP1A","SUPT6H","ZNF207","LIG3","SMARCE1","VPS25","PSME3","KPNA2","JPT1","HDGFL2",  
                                                    "RAVER1","SMARCA4","NFIX","SAMD1","UBA2","BAX","PCNA","XRN2","RBM12","RPRD1B","COL18A1","COL6A2","PRDX4","USP11","GPKOW","TCEAL4","HMGB3","BGN","HCFC1" )] <- "UP - PROTEIN"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c("ATP1A1","PBXIP1","TNR","SLC1A4","SFXN5","TUBA4A","PLCD1","LSAMP","RAB6B","PCCB","MRAS","NCEH1","PEX5L","MAP6D1","ATP5ME","SLC4A4","ENPP6","SLC25A4","TPPP","CPLANE1",  
                                                   "PURA","KCTD16","CUTA","CPNE5","RALA","GNAI1","ATP5MF","AGBL3","KBTBD11","NEFM","NEFL","ATP6V1H","SYBU","GLIPR2","NIPSNAP3A","PTGDS","CAMK2G","GLUD1","INA","CNTN1",    
                                                   "ATP5F1B","ALDH2","PEBP1","ALDH6A1","NRXN3","TEDC1","CDIP1","ABAT","GNAO1","PLLP","WDR59","CDH13","ALDOC","CNP","CNTNAP1","MBP","TUBB4A","MAG","ATP1A3","SNTA1",    
                                                   "GNAZ","NEFH","PLP1","RAP2C","MT-ATP6"  )] <- "DOWN - PROTEIN"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c("RCC1","COL6A3","RFC4","TGFBI","H1-5","MCM3", "SERPINH1","TMPO","KPNA2","JPT1","PCNA","COL18A1", "COL6A2","PRDX4","BGN")] <- "UP - Protein and RNA"
diffprotdeseq$group[diffprotdeseq$gene_name %in% c("TNR","SLC1A4","SLC4A4","ALDH2","ALDH6A1","ALDOC", "SNTA1")] <- "DOWN - Protein and RNA"



diffprotdeseq$group[diffprotdeseq$gene_name %in% c("GSTM3","ABAT","ALDH1A1","ANXA2","TNR","GOT1","KCTD21","CTNND2","HSPB1","VIM")] <- "High Correlation" 

ggplot(diffprotdeseq, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq, group == "UP - RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "UP - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "UP - Protein and RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - Protein and RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "High Correlation"), color = "black", size = 2)+
  geom_text_repel(data = subset(diffprotdeseq, group == "High Correlation"), color = "black")+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")

ggplot(diffprotdeseq, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq, group == "UP - RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "UP - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "UP - Protein and RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "DOWN - Protein and RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group == "High Correlation"), color = "black", size = 2)+
  geom_text_repel(data = subset(diffprotdeseq, group == "UP - PROTEIN"), color = "black")+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")

#HIGH COR AND UP RNA        HSPB1, VIM, ANXA2 

#HIGH COR AND DOWN RNA -    GSTM3 

#HIGH COR AND DOWN PROT -   ABAT

#HIGH COR DOWN PROTRNA -    TNR, ALDH6A1

#NOT DE -                   CTNND2, KCTD21, GOT

diffprotdeseq$group2<- "Neutral"
diffprotdeseq$group2[diffprotdeseq$gene_name %in% c ("PCNA","MCM7","NFIX","HMGN2","KPNA2","TGFBI","MDC1","RPA3","ASNS","SMARCD1")] <- "UP - PROTEIN"
diffprotdeseq$group2[diffprotdeseq$gene_name %in% c( "PLP1","PLLP","MBP","MAG","TNR","NEFL","PTGDS","WDR59","ENPP6","CNP")] <- "DOWN - PROTEIN"
diffprotdeseq$group2[diffprotdeseq$gene_name %in% c("RCC1","COL6A3","RFC4","TGFBI","H1-5","MCM3", "SERPINH1","TMPO","KPNA2","JPT1","PCNA","COL18A1", "COL6A2","PRDX4","BGN")] <- "UP - Protein and RNA"
diffprotdeseq$group2[diffprotdeseq$gene_name %in% c("TNR","SLC1A4","SLC4A4","ALDH2","ALDH6A1","ALDOC", "SNTA1")] <- "DOWN - Protein and RNA"


ggplot(diffprotdeseq, aes(x = stat, y = t, color = group2, label = gene_name))+
  geom_point(data = subset(diffprotdeseq, group2 == "Neutral"), size=0.5)+
  geom_point(data = subset(diffprotdeseq, group2 == "UP - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group2 == "DOWN - PROTEIN"), size=2)+
  geom_point(data = subset(diffprotdeseq, group2 == "UP - Protein and RNA"), size=2)+
  geom_point(data = subset(diffprotdeseq, group2 == "DOWN - Protein and RNA"), size=2)+
  geom_text_repel(data = subset(diffprotdeseq, group2 == "UP - PROTEIN"), color = "black")+
  geom_text_repel(data = subset(diffprotdeseq, group2 == "DOWN - PROTEIN"), color = "black")+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")


#Kijken of het met RNA reads te maken heeft
readcounts2<- readcount.vst %>%
  dplyr::mutate(sumrow = rowSums(.))  %>% 
  rownames_to_column("gene") 

library("writexl")
write_xlsx(diffprot2 , "~diffprot2.xlsx")

write_xlsx(deseq2Results , "~deseq2Results")

upprotRNA <- diffprot2signup %>% 
  inner_join(signupRNA, by = c("Protein" = "gene_name"))

downprotRNA <- diffprot2signdown %>% 
  inner_join(signdownRNA, by = c("Protein" = "gene_name"))


######################################PLOTSIGNUP RNA PROT
signupprotRNAproteinplot <- left_join(pre_protein_counts_filtered2, metadata_prot_RNA2, by = c("patient_id" = "Sample_Name")) %>% 
  dplyr::select("Customer_ID","PCNA", "KPNA2", "TGFBI","H1-5", "COL6A2", "COL18A1", "SERPINH1","MCM3","PRDX4","RFC4","TMPO","COL6A3", "JPT1", "RCC1", "BGN") %>% 
  dplyr::left_join(metadataALL %>% dplyr::select(Sample_Name, methylation.sub.diagnosis), by = c("Customer_ID" = "Sample_Name")) %>% 
  dplyr::filter(methylation.sub.diagnosis == "A_IDH" | methylation.sub.diagnosis =="A_IDH_HG") %>% 
  group_by(methylation.sub.diagnosis)

signupprotRNAproteinplot$methylation.sub.diagnosis <- factor(signupprotRNAproteinplot$methylation.sub.diagnosis, levels = c("A_IDH","A_IDH_HG"))

ggboxplot(signupprotRNAproteinplot, x = "methylation.sub.diagnosis", y = "PCNA", add = "dotplot", color = "methylation.sub.diagnosis", palette =c("#00AFBB", "#E7B800"))

library(ggbeeswarm)
upprotRNA$Protein[upprotRNA$Protein == "H1-5"] <-"H1_5"

signupprotRNAproteinplot<- signupprotRNAproteinplot %>% 
  dplyr::rename("H1_5" = "H1-5")

plotlistupALL = list()
for(i in 1:length(upprotRNA$Protein)){
  p <- ggboxplot(signupprotRNAproteinplot, x = "methylation.sub.diagnosis", y = upprotRNA$Protein[i], title = upprotRNA$Protein[i],
                 color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
  p<- p +  geom_quasirandom()
  
  plotlistupALL[[i]]<- p
}

ggexport (
  plotlist = plotlistupALL, filename = "upALL.pdf",
  ncol = 3, nrow = 5
)

upprotRNA$Protein[upprotRNA$Protein == "H1_5"] <-"H1-5"
signupprotRNAproteinplot<- signupprotRNAproteinplot %>% 
  dplyr::rename("H1-5" = "H1_5")

rm(plotlistupALL, signupprotRNAproteinplot)

signdownprotRNAproteinplot <- left_join(pre_protein_counts_filtered2, metadata_prot_RNA2, by = c("patient_id" = "Sample_Name")) %>% 
  dplyr::select("Customer_ID","TNR", "SLC1A4", "SLC4A4","ALDH2", "ALDH6A1", "ALDOC", "SNTA1") %>% 
  dplyr::left_join(metadataALL %>% dplyr::select(Sample_Name, methylation.sub.diagnosis), by = c("Customer_ID" = "Sample_Name")) %>% 
  dplyr::filter(methylation.sub.diagnosis == "A_IDH" | methylation.sub.diagnosis =="A_IDH_HG") %>% 
  group_by(methylation.sub.diagnosis)

signdownprotRNAproteinplot$methylation.sub.diagnosis <- factor(signdownprotRNAproteinplot$methylation.sub.diagnosis, levels = c("A_IDH","A_IDH_HG"))

plotlistdownALL = list()
for(i in 1:length(downprotRNA$Protein)){
  p <- ggboxplot(signdownprotRNAproteinplot, x = "methylation.sub.diagnosis", y = downprotRNA$Protein[i], title = downprotRNA$Protein[i],
                 color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
  p<- p +  geom_quasirandom()
  
  plotlistdownALL[[i]]<- p
}


ggexport (
  plotlist = plotlistdownALL, filename = "downALL.pdf",
  ncol = 3, nrow = 3
)

rm(plotlistdownALL, signdownprotRNAproteinplot)


inputplotRNAdiffexpr <- readcount.vst %>% 
  rownames_to_column("gene_id") %>% 
  dplyr::filter(gene_id %in% c("ENSG00000132646.11", "ENSG00000182481.9", "ENSG00000120708.17", "ENSG00000184357.5",
                "ENSG00000142173.16", "ENSG00000182871.16", "ENSG00000149257.15", "ENSG00000112118.20",
                "ENSG00000123131.13", "ENSG00000163918.10", "ENSG00000120802.13", "ENSG00000163359.16",
                "ENSG00000189159.16", "ENSG00000180198.16", "ENSG00000182492.16")) %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate(gene_id = NULL) %>% 
  select(gene_name, everything()) %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  left_join(metadata %>% dplyr::select(GS_ID, methylation.sub.diagnosis))

inputplotRNAdiffexpr$methylation.sub.diagnosis <- factor(inputplotRNAdiffexpr$methylation.sub.diagnosis, levels = c("A_IDH","A_IDH_HG"))


upprotRNA$Protein[upprotRNA$Protein == "H1-5"] <-"H1_5"

inputplotRNAdiffexpr<- inputplotRNAdiffexpr %>% 
  dplyr::rename("H1_5" = "H1-5")

plotlistupALLRNA = list()
for(i in 1:length(upprotRNA$Protein)){
  p <- ggboxplot(inputplotRNAdiffexpr, x = "methylation.sub.diagnosis", y = upprotRNA$Protein[i], title = upprotRNA$Protein[i],
                 color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
  p<- p +  geom_quasirandom()
  
  plotlistupALLRNA[[i]]<- p
}

ggexport (
  plotlist = plotlistupALLRNA, filename = "upALLRNA.pdf",
  ncol = 3, nrow = 5
)

upprotRNA$Protein[upprotRNA$Protein == "H1_5"] <-"H1-5"

rm(inputplotRNAdiffexpr, plotlistupALLRNA)

inputplotRNAdiffexprdown <- readcount.vst %>% 
  rownames_to_column("gene_id") %>% 
  dplyr::filter(gene_id %in% c("ENSG00000116147.17", "ENSG00000115902.11", "ENSG00000111275.13", "ENSG00000101400.6",
                               "ENSG00000080493.17", "ENSG00000109107.14", "ENSG00000119711.13")) %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate(gene_id = NULL) %>% 
  select(gene_name, everything()) %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID") %>% 
  left_join(metadata %>% dplyr::select(GS_ID, methylation.sub.diagnosis))

inputplotRNAdiffexprdown$methylation.sub.diagnosis <- factor(inputplotRNAdiffexprdown$methylation.sub.diagnosis, levels = c("A_IDH","A_IDH_HG"))

plotlistdownRNAALL = list()
for(i in 1:length(downprotRNA$Protein)){
  p <- ggboxplot(inputplotRNAdiffexprdown, x = "methylation.sub.diagnosis", y = downprotRNA$Protein[i], title = downprotRNA$Protein[i],
                 color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
  p<- p +  geom_quasirandom()
  
  plotlistdownRNAALL[[i]]<- p
}


ggexport (
  plotlist = plotlistdownRNAALL, filename = "downRNAALL.pdf",
  ncol = 3, nrow = 3
)