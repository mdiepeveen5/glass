library("ggplot2")
library("ggrepel")
#needs diffprot2 and deseq2Results and gene.annot2

load(file = "diffprot2.Rdata")
load(file = "deseq2Results.Rdata")

diffprotdeseq2 <- deseq2Results %>% 
  left_join(gene.annot2) %>% 
  inner_join(diffprot2, by = c("gene_name"="Protein"))


diffprotdeseq2$group<- "Not differentially expressed"


diffprotdeseq2$group[diffprotdeseq2$gene_name %in% DNAbindingdeseq2$gene] <- "DNA binding"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% cycling.cell.markers$gene_name] <- "Cell cycle"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("H4C1", "H3-3A", "H2BC17", "H2BC12", "H2AZ1","H2AC4", "H2AC21", "H2AC20", "H2AC11", "H1-5", "H1-4", "H1-3","H1-2", "H1-10","H1-0")] <- "Histones"


#DONT USE THIS
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c ("RCC2","SRRM1","HMGN2","RCC1","RBBP4","ANP32E","EEF1AKNMT","CDC73","H3-3A","PRPF40A","CNOT9","COL6A3","ACTL6A","SOX2","RFC4","PDS5A","FIP1L1","PPAT","SMARCA5","TGFBI",   
                                                      "TXNDC5","DEK","H1-5","MDC1","MAPK14","MCM3","RPF2","NUP43","RPA3","CBX3","ASNS","BUD31","MCM7","ENY2","RAD21","NCBP1","PRPF4","GOLGA2","GTF3C5","DDX50",   
                                                      "SSRP1","MTA2","SERPINH1","MRE11","DDX47","SMARCD1","SNRPF","TMPO","NOVA1","RSL1D1","MYBBP1A","SUPT6H","ZNF207","LIG3","SMARCE1","VPS25","PSME3","KPNA2","JPT1","HDGFL2",  
                                                      "RAVER1","SMARCA4","NFIX","SAMD1","UBA2","BAX","PCNA","XRN2","RBM12","RPRD1B","COL18A1","COL6A2","PRDX4","USP11","GPKOW","TCEAL4","HMGB3","BGN","HCFC1" )] <- "Differential up"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c ("PCNA", "MCM7", "MDC1", "RPA3", "ASNS", "RPRD1B","CDC73", "SOX2", "MYBBP1A",
                                                      "SMARCA5", "GOLGA2", "RAD21", "MCM3", "ZNF207", "PDS5A", "PSME3", "ACTL6A",
                                                      "MAPK14", "RCC2","LIG3","MRE11", "RBBP4", "RCC1","NUP43","PRPF40A","HCFC1",
                                                      "BAX", "EEF1AKNMT")] <- "Cell Cycle"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("ATP1A1","PBXIP1","TNR","SLC1A4","SFXN5","TUBA4A","PLCD1","LSAMP","RAB6B","PCCB","MRAS","NCEH1","PEX5L","MAP6D1","ATP5ME","SLC4A4","ENPP6","SLC25A4","TPPP","CPLANE1",  
                                                     "PURA","KCTD16","CUTA","CPNE5","RALA","GNAI1","ATP5MF","AGBL3","KBTBD11","NEFM","NEFL","ATP6V1H","SYBU","GLIPR2","NIPSNAP3A","PTGDS","CAMK2G","GLUD1","INA","CNTN1",    
                                                     "ATP5F1B","ALDH2","PEBP1","ALDH6A1","NRXN3","TEDC1","CDIP1","ABAT","GNAO1","PLLP","WDR59","CDH13","ALDOC","CNP","CNTNAP1","MBP","TUBB4A","MAG","ATP1A3","SNTA1",    
                                                     "GNAZ","NEFH","PLP1","RAP2C","MT-ATP6"  )] <- "Differential down"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("PCNA", "MCM3","MCM7", "RCC1")] <- "Paralogs upgerulated in RNA"
diffprotdeseq2$group[diffprotdeseq2$gene_name %in% c("H4C1", "H3-3A", "H2BC17", "H2BC12", "H2AZ1","H2AC4", "H2AC21", "H2AC20", "H2AC11", "H1-5", "H1-4", "H1-3","H1-2", "H1-10","H1-0")] <- "Histones"




ggplot(diffprotdeseq2, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq2, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq2, group == "Differential up"), size=1)+
  geom_point(data = subset(diffprotdeseq2, group == "Differential down"), size=1)+
  geom_point(data = subset(diffprotdeseq2, group == "Histones"), size=3)+
  geom_point(data = subset(diffprotdeseq2, group == "Cell cycle"), size=3)+
  geom_point(data = subset(diffprotdeseq2, group == "Paralogs upgerulated in RNA"), size=3)+
  geom_text_repel(data = subset(diffprotdeseq2, group == "Cell cycle"), color = "black")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_abline(intercept = 0, slope = 0.5)+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")+
  scale_color_manual(values = c("grey", "blue", "red", "brown", "mediumpurple", "yellow"),
                     name = "Group",
                     breaks=c("Not differentially expressed", "Differential up", "Differential down", "Histones", "Cell Cycle",
                              "Paralogs upgerulated in RNA"),
                     labels = c("Not differentially expressed", "Differential up", "Differential down", "Histones", "Cell Cycle",
                                "Paralogs upgerulated in RNA"))



ggplot(diffprotdeseq2, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq2, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq2, group == "DNA binding"), size=2)+
  geom_point(data = subset(diffprotdeseq2, group == "Histones"), size=3)+
  geom_point(data = subset(diffprotdeseq2, group == "Cell cycle"), size=3)+
  geom_text_repel(data = subset(diffprotdeseq2, group == "Cell cycle"), color = "black")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_abline(intercept = 0, slope = 0.5)+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")+
  scale_color_manual(values = c("grey", "royalblue", "brown", "purple"),
                     name = "Group",
                     breaks=c("Not differentially expressed", "DNA binding", "Histones", "Cell cycle"),
                     labels = c("Not differentially expressed", "DNA binding", "Histones", "Cell cycle"))


#t-test proteomics
t.test(subset(diffprotdeseq2, group == "Cell cycle") %>% dplyr::select(stat), subset(diffprotdeseq2, group != "Cell cycle") %>% dplyr::select(stat), var.equal = T)

#p-value < 2.2e-16

t.test(subset(diffprotdeseq2, group == "DNA binding") %>% dplyr::select(stat), subset(diffprotdeseq2, group != "DNA binding") %>% dplyr::select(stat), var.equal = T)

#p-value = 2.535e-07

#t-test RNA

t.test(subset(diffprotdeseq2, group == "Cell cycle") %>% dplyr::select(t), subset(diffprotdeseq2, group != "Cell cycle") %>% dplyr::select(t), var.equal = T)

#p-vale = 4.251e-07

t.test(subset(diffprotdeseq2, group == "DNA binding") %>% dplyr::select(t), subset(diffprotdeseq2, group != "DNA binding") %>% dplyr::select(t), var.equal = T)

#p-value < 2.2e-16




##GRADE


diffprotdeseq3 <- deseq2ResultsGRADE %>% 
  left_join(gene.annot2) %>% 
  inner_join(diffprot3, by = c("gene_name"="Protein"))


diffprotdeseq3$group<- "Not differentially expressed"


diffprotdeseq3$group[diffprotdeseq3$gene_name %in% DNAbindingdeseq2$gene] <- "DNA binding"
diffprotdeseq3$group[diffprotdeseq3$gene_name %in% cycling.cell.markers$gene_name] <- "Cell cycle"
diffprotdeseq3$group[diffprotdeseq3$gene_name %in% c("H4C1", "H3-3A", "H2BC17", "H2BC12", "H2AZ1","H2AC4", "H2AC21", "H2AC20", "H2AC11", "H1-5", "H1-4", "H1-3","H1-2", "H1-10","H1-0")] <- "Histones"


ggplot(diffprotdeseq3, aes(x = stat, y = t, color = group, label = gene_name))+
  geom_point(data = subset(diffprotdeseq3, group == "Not differentially expressed"), size=0.5)+
  geom_point(data = subset(diffprotdeseq3, group == "DNA binding"), size=2)+
  geom_point(data = subset(diffprotdeseq3, group == "Histones"), size=3)+
  geom_point(data = subset(diffprotdeseq3, group == "Cell cycle"), size=3)+
  geom_text_repel(data = subset(diffprotdeseq3, group == "Cell cycle"), color = "black")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_abline(intercept = 0, slope = 0.5)+
  theme_bw()+
  xlab("RNA")+
  ylab("Protein")+
  scale_color_manual(values = c("grey", "royalblue", "brown", "purple"),
                     name = "Group",
                     breaks=c("Not differentially expressed", "DNA binding", "Histones", "Cell cycle"),
                     labels = c("Not differentially expressed", "DNA binding", "Histones", "Cell cycle"))
