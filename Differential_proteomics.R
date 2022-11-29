# #neem load_proteomics_data
# diffproteomics <- read.csv("~/projects/glass/Limma_DEP_result.csv") %>% 
#   dplyr::rename(Protein = X)
# 
# res.paired.a.exon <-readRDS("~/projects/glass/res.paired.a.exon.Rds")
# 
# diffproteomics$Protein[diffproteomics$Protein == "WARS"] <-"WARS1"
# diffproteomics$Protein[diffproteomics$Protein == "NARS"] <-"NARS1"
# diffproteomics$Protein[diffproteomics$Protein == "HARS"] <-"HARS1"
# diffproteomics$Protein[diffproteomics$Protein == "DARS"] <-"DARS1"
# diffproteomics$Protein[diffproteomics$Protein == "TARS"] <-"TARS1"
# diffproteomics$Protein[diffproteomics$Protein == "VARS"] <-"VARS1"
# diffproteomics$Protein[diffproteomics$Protein == "GARS"] <-"GARS1"
# diffproteomics$Protein[diffproteomics$Protein == "IARS"] <-"IARS1"
# diffproteomics$Protein[diffproteomics$Protein == "QARS"] <-"QARS1"
# diffproteomics$Protein[diffproteomics$Protein == "AARS"] <-"AARS1"
# diffproteomics$Protein[diffproteomics$Protein == "CARS"] <-"CARS1"
# diffproteomics$Protein[diffproteomics$Protein == "SARS"] <-"SARS1"
# diffproteomics$Protein[diffproteomics$Protein == "RARS"] <-"RARS1"
# diffproteomics$Protein[diffproteomics$Protein == "YARS"] <-"YARS1"
# diffproteomics$Protein[diffproteomics$Protein == "MARS"] <-"MARS1"
# diffproteomics$Protein[diffproteomics$Protein == "KARS"] <-"KARS1"
# diffproteomics$Protein[diffproteomics$Protein == "LARS"] <-"LARS1"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H2BK"] <-"H2BC12"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H2AB"] <-"H2AC4"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H2AG"] <-"H2AC11"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H1E"] <-"H1-4"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H1B"] <-"H1-5"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H1D"] <-"H1-3"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H1C"] <-"H1-2"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H2BO"] <-"H2BC17"
# diffproteomics$Protein[diffproteomics$Protein == "HIST1H4A"] <-"H4C1"
# diffproteomics$Protein[diffproteomics$Protein == "H3F3A"] <-"H3-3A"
# diffproteomics$Protein[diffproteomics$Protein == "HIST2H2AC"] <-"H2AC20"
# diffproteomics$Protein[diffproteomics$Protein == "HIST2H2AB"] <-"H2AC21"
# diffproteomics$Protein[diffproteomics$Protein == "H1FX"] <-"H1-10"
# diffproteomics$Protein[diffproteomics$Protein == "H2AFY2"] <-"MACROH2A2"
# diffproteomics$Protein[diffproteomics$Protein == "H1F0"] <-"H1-0"
# diffproteomics$Protein[diffproteomics$Protein == "H2AFY"] <-"MACROH2A1"
# diffproteomics$Protein[diffproteomics$Protein == "H2AFZ"] <-"H2AZ1"
# 
# 
# 
# up1 = dge.partially.paired.clusters %>% 
#   dplyr::filter(up.1 == T) %>% 
#   dplyr::select(gene_name) %>% 
#   dplyr::inner_join(interprotRNA)
# 
# up2 = dge.partially.paired.clusters %>% 
#   dplyr::filter(up.2 == T) %>% 
#   dplyr::select(gene_name) %>% 
#   dplyr::inner_join(interprotRNA)
# 
# up3 = dge.partially.paired.clusters %>% 
#   dplyr::filter(up.3 == T) %>% 
#   dplyr::select(gene_name) %>% 
#   dplyr::inner_join(interprotRNA)
# 
# #UP1
# up1diffprot <- corproteinorderedNATR %>% 
#   dplyr::select("LMNB1", "ACTC1", "CDK1", "H4C1", "H2AC4", "H2AC11", "H2BC12", "H1-5", "H2BC17") %>% 
#   tibble::rownames_to_column("GS_ID") %>% 
#   dplyr::left_join(metadata_prot_RNA %>% dplyr::select(GS_ID, resec2), by = c("GS_ID")) %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up1diffprot$colors <- factor(up1diffprot$colors, levels = c("Primary","Recurrent"))
# 
# plotlistup1 = list()
# for(i in 1:length(up1$gene_name)){
#   p <- ggboxplot(up1diffprot, x = "colors", y = up1$gene_name[i], title = up1$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
#   plotlistup1[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup1, filename = "up1.pdf",
#   ncol = 3, nrow = 3
# )
# 
# rm(plotlistup1, p)
# 
# #ggboxplot(up1diffprot, x= "colors", y = "LMNB1",
# #          color = "colors", palette = c("#00AFBB", "#E7B800"),
# #          title = "LMNB1",
# #          ylab = "log2(Intensity)", xlab = "Groups",)
# 
# #UP2
# 
# up2diffprot <- corproteinorderedNATR %>% 
#   dplyr::select("HSPG2", "CHI3L1", "COL6A3","TGFBI", "COL1A2", "VGF", "LUM","COL4A2","ACAN","COL1A1","FCGBP","COL6A2") %>% 
#   tibble::rownames_to_column("GS_ID") %>% 
#   dplyr::left_join(metadata_prot_RNA %>% dplyr::select(GS_ID, resec2), by = c("GS_ID")) %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up2diffprot$colors <- factor(up2diffprot$colors, levels = c("Primary","Recurrent"))
# 
# plotlistup2 = list()
# for(i in 1:length(up2$gene_name)){
#   p <- ggboxplot(up2diffprot, x = "colors", y = up2$gene_name[i], title = up2$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
#   plotlistup2[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup2, filename = "up2.pdf",
#   ncol = 3, nrow = 4
# )
# 
# rm(plotlistup2, p)
# 
# #UP3
# up3diffprot <- corproteinorderedNATR %>% 
#   dplyr::select("AQP1", "FABP5") %>% 
#   tibble::rownames_to_column("GS_ID") %>% 
#   dplyr::left_join(metadata_prot_RNA %>% dplyr::select(GS_ID, resec2), by = c("GS_ID")) %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up3diffprot$colors <- factor(up3diffprot$colors, levels = c("Primary","Recurrent"))
# 
# plotlistup3 = list()
# for(i in 1:length(up3$gene_name)){
#   p <- ggboxplot(up3diffprot, x = "colors", y = up3$gene_name[i], title = up3$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
#   plotlistup3[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup3, filename = "up3.pdf",
#   ncol = 2, nrow = 1
# )
# 
# rm(plotlistup3, p)
# 
# 
# primary <- dplyr::filter(up1diffprot, colors == "Primary" )
# recurrent <- dplyr::filter(up1diffprot, colors=="Recurrent")
# 
# cor.test(primary$LMNB1, recurrent$LMB1)
# 
# cor.test(dplyr::filter(up1diffprot, colors == "Primary" ) %>% dplyr::filter(up1diffprot, colors == "Recurrent" ) %>% dplyr::select("LMNB1"))
# 

##NOT FILTERED PROTEIN DATA
tmpprot <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) %>%
  dplyr::rename(Protein = X)
tmpprot$Protein[tmpprot$Protein == "WARS"] <-"WARS1"
tmpprot$Protein[tmpprot$Protein == "NARS"] <-"NARS1"
tmpprot$Protein[tmpprot$Protein == "HARS"] <-"HARS1"
tmpprot$Protein[tmpprot$Protein == "DARS"] <-"DARS1"
tmpprot$Protein[tmpprot$Protein == "TARS"] <-"TARS1"
tmpprot$Protein[tmpprot$Protein == "VARS"] <-"VARS1"
tmpprot$Protein[tmpprot$Protein == "GARS"] <-"GARS1"
tmpprot$Protein[tmpprot$Protein == "IARS"] <-"IARS1"
tmpprot$Protein[tmpprot$Protein == "QARS"] <-"QARS1"
tmpprot$Protein[tmpprot$Protein == "AARS"] <-"AARS1"
tmpprot$Protein[tmpprot$Protein == "CARS"] <-"CARS1"
tmpprot$Protein[tmpprot$Protein == "SARS"] <-"SARS1"
tmpprot$Protein[tmpprot$Protein == "RARS"] <-"RARS1"
tmpprot$Protein[tmpprot$Protein == "YARS"] <-"YARS1"
tmpprot$Protein[tmpprot$Protein == "MARS"] <-"MARS1"
tmpprot$Protein[tmpprot$Protein == "KARS"] <-"KARS1"
tmpprot$Protein[tmpprot$Protein == "LARS"] <-"LARS1"
tmpprot$Protein[tmpprot$Protein == "HIST1H2BK"] <-"H2BC12"
tmpprot$Protein[tmpprot$Protein == "HIST1H2AB"] <-"H2AC4"
tmpprot$Protein[tmpprot$Protein == "HIST1H2AG"] <-"H2AC11"
tmpprot$Protein[tmpprot$Protein == "HIST1H1E"] <-"H1-4"
tmpprot$Protein[tmpprot$Protein == "HIST1H1B"] <-"H1-5"
tmpprot$Protein[tmpprot$Protein == "HIST1H1D"] <-"H1-3"
tmpprot$Protein[tmpprot$Protein == "HIST1H1C"] <-"H1-2"
tmpprot$Protein[tmpprot$Protein == "HIST1H2BO"] <-"H2BC17"
tmpprot$Protein[tmpprot$Protein == "HIST1H4A"] <-"H4C1"
tmpprot$Protein[tmpprot$Protein == "H3F3A"] <-"H3-3A"
tmpprot$Protein[tmpprot$Protein == "HIST2H2AC"] <-"H2AC20"
tmpprot$Protein[tmpprot$Protein == "HIST2H2AB"] <-"H2AC21"
tmpprot$Protein[tmpprot$Protein == "H1FX"] <-"H1-10"
tmpprot$Protein[tmpprot$Protein == "H2AFY2"] <-"MACROH2A2"
tmpprot$Protein[tmpprot$Protein == "H1F0"] <-"H1-0"
tmpprot$Protein[tmpprot$Protein == "H2AFY"] <-"MACROH2A1"
tmpprot$Protein[tmpprot$Protein == "H2AFZ"] <-"H2AZ1"
tmpprot$Protein[tmpprot$Protein == "TARSL2"] <-"TARS3"
tmpprot$Protein[tmpprot$Protein == "ASNA1"] <-"GET3"
tmpprot$Protein[tmpprot$Protein == "EPRS"] <-"EPRS1"
tmpprot$Protein[tmpprot$Protein == "ADSS"] <-"ADSS2"
tmpprot$Protein[tmpprot$Protein == "CRAD"] <-"CRACD"
tmpprot$Protein[tmpprot$Protein == "HIST1H2BA"] <-"H2BC1"
tmpprot$Protein[tmpprot$Protein == "KIF1BP"] <-"KIFBP"
tmpprot$Protein[tmpprot$Protein == "FAM129B"] <-"NIBAN2"
tmpprot$Protein[tmpprot$Protein == "FAM49A"] <-"CYRIA"
tmpprot$Protein[tmpprot$Protein == "RYDEN"] <-"SHFL"
tmpprot$Protein[tmpprot$Protein == "FAM49B"] <-"CYRIB"
tmpprot$Protein[tmpprot$Protein == "ADPRHL2"] <-"ADPRS"

pre_protein_counts_filtered2 <- tmpprot %>%
  tibble::column_to_rownames("Protein") %>%
  t() %>%
  as.data.frame()%>%
  tibble::rownames_to_column("patient_id")


metadata_prot_RNA2 <- data.frame(Sample_Name = colnames(tmpprot[,-1])) %>%
  dplyr::mutate(sid = gsub('^[^_]+_[^_]+_([0-9]+)_.+$', '\\1', Sample_Name)) %>%
  dplyr::mutate(resec = gsub('^[^_]+_[^_]+_[0-9]+_(.+)$', '\\1', Sample_Name)) %>%
  dplyr::mutate(resec2 = case_when(
    resec == "P" ~ "R1",
    resec == "R1" ~ "R2",
    resec == "R2" ~ "R2",
    TRUE ~ "#$@#$"
  )) %>%
  dplyr::mutate(Customer_ID=paste0(sid, "_", resec2)) %>%
  dplyr::mutate(resec = NULL) %>%
  dplyr::mutate(sid = NULL)
# 
# 
# rm(tmpprot)
# 
# up1diffprotALL <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id", "LMNB1", "H4C1", "H2AC4", "H2AC11", "H2BC12","H1-5", "H2BC17","CDK1", "ACTC1") %>% 
#   dplyr::left_join(metadata_prot_RNA2 %>% dplyr::select(Sample_Name, resec2), by = c("patient_id" = "Sample_Name")) %>% 
#   dplyr::filter(resec2 != "R3") %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up1diffprotALL$colors <- factor(up1diffprotALL$colors, levels = c("Primary","Recurrent"))
# 
# up1$gene_name[up1$gene_name == "H1-5"] <-"H1_5"
# 
# up1diffprotALL<- up1diffprotALL %>% 
#   dplyr::rename("H1_5" = "H1-5")
# 
# plotlistup1ALL = list()
# for(i in 1:length(up1$gene_name)){
#   p <- ggboxplot(up1diffprotALL, x = "colors", y = up1$gene_name[i], title = up1$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups")
#   
#   p<- p +  geom_quasirandom()
#   
#   plotlistup1ALL[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup1ALL, filename = "up1ALL.pdf",
#   ncol = 3, nrow = 3
# )
# rm(plotlistup1ALL, up1diffprotALL)
# 
# up1$gene_name[up1$gene_name == "H1_5"] <-"H1-5"
# 
# ggboxplot(up1diffprotALL, x = "colors", y = "H2AC11", add = "dotplot")
# 
# #UP2
# up2diffprotALL <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id","HSPG2", "CHI3L1", "COL6A3","TGFBI", "COL1A2", "VGF", "LUM","COL4A2","ACAN","COL1A1","FCGBP","COL6A2") %>% 
#   dplyr::left_join(metadata_prot_RNA2 %>% dplyr::select(Sample_Name, resec2), by = c("patient_id" = "Sample_Name")) %>% 
#   dplyr::filter(resec2 != "R3") %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up2diffprotALL$colors <- factor(up2diffprotALL$colors, levels = c("Primary","Recurrent"))
# 
# plotlistup2ALL = list()
# for(i in 1:length(up2$gene_name)){
#   p <- ggboxplot(up2diffprotALL, x = "colors", y = up2$gene_name[i], title = up2$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups")
#   p<- p +  geom_quasirandom()
#   plotlistup2ALL[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup2ALL, filename = "up2ALL.pdf",
#   ncol = 3, nrow = 4
# )
# 
# rm(up2diffprotALL)
# 
# 
# #UP3
# 
# up3diffprotALL <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id", "AQP1", "FABP5") %>% 
#   dplyr::left_join(metadata_prot_RNA2 %>% dplyr::select(Sample_Name, resec2), by = c("patient_id" = "Sample_Name")) %>% 
#   dplyr::mutate(colors = ifelse(resec2 == "R1", "Primary",  "Recurrent")) %>% 
#   group_by(colors)
# 
# up3diffprotALL$colors <- factor(up3diffprotALL$colors, levels = c("Primary","Recurrent"))
# 
# plotlistup3ALL = list()
# for(i in 1:length(up3$gene_name)){
#   p <- ggboxplot(up3diffprotALL, x = "colors", y = up3$gene_name[i], title = up3$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups",)
#   p<- p +  geom_quasirandom()
#   plotlistup3ALL[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup3ALL, filename = "up3ALL.pdf",
#   ncol = 2, nrow = 1
# )
# 
# rm(up3diffprotALL, p)
# 
# 
# #LIMMA
# 
# library(limma)
# 
diffprotWHO <- diffprotWHO %>% 
   left_join(metadata_prot_RNA2)

DNA_binding <- read_csv("DNA_binding.csv")
 
 table(metadata_prot_RNA2$resec2)
 metadata_prot_RNA2$resec2
 
 y2 <- pre_protein_counts_filtered2 %>% 
   tibble::column_to_rownames("patient_id") %>% 
   t() %>% 
   as.data.frame()
 
 design2 <- model.matrix(~metadata_prot_RNA2$resec2)
 
 
 fit <- lmFit(y2, design2)
 fit <- eBayes(fit)
 diffprot <- topTable(fit,sort="none",n=Inf) 
 
 diffprot <- diffprot %>% 
   rownames_to_column("Protein")
 
 volcanodiffprotprimrecur <- diffprot %>% 
   rownames_to_column("gene") %>% 
   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75,  "Significant", "Not Significant")) %>% 
   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & gene %in% DNA_binding$name, "DNA binding", colors)) %>% 
   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & gene %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>% 
   dplyr::mutate(size = ifelse(colors == "DNA binding" | colors == "Cell cycle", 1.5,1))
 
 
pdf("Volcano_prot_primrecur.pdf", width = 8, height = 5)
 
ggplot(volcanodiffprotprimrecur, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
   geom_point(aes(col=colors, size = size)) + 
   scale_radius(range = c(1,2))+
   guides(size = "none")+
   scale_color_manual(values=c("grey")) + 
   theme_bw()+
   geom_text_repel(data=subset(volcanodiffprotprimrecur, colors == "Cell cycle" | colors == "DNA binding"), aes(label=gene)) +
   geom_hline(yintercept=1.3, linetype="dashed")+
   geom_vline(xintercept = 0.75, linetype="dashed")+
   geom_vline(xintercept = -0.75, linetype="dashed")+
   labs(title = 'Differentially expressed genes for proteomics data', subtitle = 'Between primary (n = 26) and recurrent (n= 29) samples')+
   xlab('Log fold change')+
   ylab('-log10(adjusted p-value)')
 
dev.off()

 rm(fit, y2, design2)
# 
# 
# 
# up3diffprotALL2 <- up3diffprotALL %>% 
#   dplyr::filter(resec2 == "R1")
# 
# mean(up3diffprotALL2$AQP1)
# 
# 
# up3diffprotALL3 <- up3diffprotALL %>% 
#   dplyr::filter(resec2 == "R2")
# mean(up3diffprotALL3$AQP1)

#######################DATA INLEZEN YOURI

library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
load("metadataALL.Rdata")
diffprotWHO <- metadata_prot_RNA2 %>% 
  dplyr::select(Customer_ID, Sample_Name) %>% 
  left_join(metadataALL, by = c("Customer_ID" = "Sample_Name"))
#  dplyr::filter(Sample_Name %in% c("AC_GII_146_P", "AC_GII_131_P", "GB_GIV_153_R1", "AC_GII_136_P", "AC_GII_141_R1") == F) #maybe needed for methylation fit but not for grade vs grade




annotations <- readxl::read_xlsx('~/projects/glass//data//glass/Proteomics/2022-03-31_data_update/Annotation_Reduced_withControls.xlsx')

annotations_filtered <- inner_join(diffprotWHO, annotations, by = c("Customer_ID"="Sample_Name"))

annotations_filtered <- annotations_filtered %>% 
  dplyr::mutate_all(function(x){return(ifelse(x == "NA",is.na(x),x))}) %>% 
  transform(OS_diagnosis_months = as.numeric(OS_diagnosis_months),
            OS_diagnosis_years = as.numeric(OS_diagnosis_years),
            OS_FirstSurgery_months = as.numeric(OS_FirstSurgery_months),
            OS_FirstSurgery_years= as.numeric(OS_FirstSurgery_years),
            Survival_FirstRecurrenceSurgery_months = as.numeric(Survival_FirstRecurrenceSurgery_months),
            Survival_FirstRecurrenceSurgery_years = as.numeric(Survival_FirstRecurrenceSurgery_years),
            PFS_years = as.numeric(PFS_years),
            PFS_months = as.numeric(PFS_months),
            Age.at.diagnosis = as.numeric(Age.at.diagnosis),
            TumorPercentage = as.numeric(TumorPercentage),
            Deceased = as.numeric(Deceased))


rm(tmpprot)

save(diffprotWHO2, file = "diffprotWHO2.Rdata")



table(diffprotWHO$methylation.sub.diagnosis)

diffprotWHO <- diffprotWHO %>%
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Customer_ID))

y3 <- pre_protein_counts_filtered2 %>% 
  tibble::column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(GB_GIV_153_R1 = NULL) %>% 
  dplyr::mutate(AC_GII_131_P = NULL) %>% 
  dplyr::mutate(AC_GII_146_P = NULL) %>% 
  dplyr::mutate(AC_GII_136_P = NULL) %>% 
  dplyr::mutate(AC_GII_141_R1 = NULL)

  
diffprotWHO2 <- diffprotWHO %>% 
  dplyr::filter(methylation.sub.diagnosis %in% c("O_IDH", "CONTR_HEMI") == F) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_131_P", "AC_GII_146_P", "GB_GIV_153_R1") ==F)
  
design3 <- model.matrix(~diffprotWHO2$methylation.sub.diagnosis)
head(design3)


fit <- lmFit(y3, design3)
fit <- eBayes(fit)
diffprot2 <- topTable(fit,sort="none",n=Inf) 

diffprot2 <- diffprot2 %>% 
  rownames_to_column("Protein")

diffprot2 <- diffprot2[order(diffprot2$adj.P.Val),]


save(diffprot2, file = "diffprot2.Rdata")

diffprot2signup <- diffprot2 %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC >0.75)

diffprot2signup <- diffprot2signup[order(-(diffprot2signup$logFC)),]

diffprot2signdown <- diffprot2 %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC < -0.75 )

diffprot2signdown <- diffprot2signdown[order((diffprot2signdown$logFC)),]

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

DNA_replication <- read_csv("DNA_repl.csv", show_col_types = FALSE)

oligo_develop <- read_csv("oligo_develop.csv", show_col_types = FALSE)

volcanodiffprot <- diffprot2 %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75,  "Significant", "Not Significant")) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & logFC > 0.75 & Protein %in% DNA_replication$name, "DNA replication", colors)) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & Protein %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & Protein %in% oligo_develop$name, "Oligodendrocyte development", colors)) %>%
  dplyr::mutate(size = ifelse(colors == "DNA replication" | colors == "Cell cycle" | colors == "Oligodendrocyte development", 1.5,1))


pdf("Volcano_prot.pdf", width = 10, height = 7)
ggplot(volcanodiffprot, aes(logFC, -log10(adj.P.Val), color = colors, label = Protein)) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(data = subset(volcanodiffprot, colors == "Not Significant"), size = 1)+
  geom_point(data = subset(volcanodiffprot, colors == "Significant"), size = 1)+
  geom_point(data = subset(volcanodiffprot, colors == "DNA replication"), size = 2)+
  geom_point(data = subset(volcanodiffprot, colors == "Cell cycle"), size = 2)+
  geom_point(data = subset(volcanodiffprot, colors == "Oligodendrocyte development"), size = 2)+
  scale_color_manual(values = c("gray87", "gray40", "royalblue", "purple", "red"),
                     name = "colors",
                     breaks=c("Not Significant", "Significant", "DNA replication", "Cell cycle", "Oligodendrocyte development"),
                     labels = c("Not Significant", "Significant", "DNA replication", "Cell cycle", "Oligodendrocyte development"))+
  theme_bw()+
  geom_text_repel(data=subset(volcanodiffprot, colors == "Cell cycle" | colors == "DNA replication" | colors == "Oligodendrocyte development"), aes(label=Protein), color = "black") +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for proteomics data', subtitle= 'Between methylation classifier LG (n =40) versus HG (n=10)')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-2.3,2.3)


dev.off()
# 
print(as.data.frame(diffprot2signup$Protein), row.names = FALSE)

print(as.data.frame(diffprot2signdown$Protein), row.names = FALSE)


print(as.data.frame(diffprot2$Protein[2997:3249],), row.names = FALSE)

#https://biit.cs.ut.ee/gplink/l/kgIxIeBOQa


diffprotWHO2 <- diffprotWHO2%>% dplyr::mutate(Resection = gsub(".*_", "", Customer_ID))

sum(diffprotWHO2$Resection == "R1")
sum(diffprotWHO2$Resection == "R2")

sum(diffprotWHO2$Resection == "R1" & diffprotWHO2$methylation.sub.diagnosis == "A_IDH_HG")

sum(diffprotWHO2$Resection == "R2" & diffprotWHO2$methylation.sub.diagnosis == "A_IDH")


rm(fit, y3, design3)


 
# 
# 
# #2 en 3 samen tegen 4
# #UP1
# up1diffprotALLWHO <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id", "LMNB1", "H4C1", "H2AC4", "H2AC11", "H2BC12","H1-5", "H2BC17","CDK1", "ACTC1")%>% 
#   dplyr::inner_join(diffprotWHO2 %>% dplyr::select(Sample_Name, methylation.sub.diagnosis), by = c("patient_id" = "Sample_Name")) %>% 
#   group_by(methylation.sub.diagnosis)
# 
# up1diffprotALLWHO$methylation.sub.diagnosis <- factor(up1diffprotALLWHO$methylation.sub.diagnosis, levels = c("A_IDH", "A_IDH_HG"))
# 
# 
# up1$gene_name[up1$gene_name == "H1-5"] <-"H1_5"
# 
# up1diffprotALLWHO<- up1diffprotALLWHO %>% 
#   dplyr::rename("H1_5" = "H1-5")
# 
# plotlistup1ALLWHO = list()
# for(i in 1:length(up1$gene_name)){
#   p <- ggboxplot(up1diffprotALLWHO, x = "methylation.sub.diagnosis", y = up1$gene_name[i], title = up1$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups")
#   
#   p<- p +  geom_quasirandom()
#   
#   plotlistup1ALLWHO[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup1ALLWHO, filename = "up1ALLWHO.pdf",
#   ncol = 3, nrow = 3
# )
# rm(plotlistup1ALLWHO, up1diffprotALLWHO)
# 
# up1$gene_name[up1$gene_name == "H1_5"] <-"H1-5"
# 
# #UP2
# up2diffprotALLWHO <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id","HSPG2", "CHI3L1", "COL6A3","TGFBI", "COL1A2", "VGF", "LUM","COL4A2","ACAN","COL1A1","FCGBP","COL6A2")%>% 
#   dplyr::inner_join(diffprotWHO2 %>% dplyr::select(Sample_Name, methylation.sub.diagnosis), by = c("patient_id" = "Sample_Name")) %>% 
#   group_by(methylation.sub.diagnosis)
# 
# up2diffprotALLWHO$methylation.sub.diagnosis <- factor(up2diffprotALLWHO$methylation.sub.diagnosis, levels = c("A_IDH", "A_IDH_HG"))
# 
# plotlistup2ALLWHO = list()
# for(i in 1:length(up2$gene_name)){
#   p <- ggboxplot(up2diffprotALLWHO, x = "methylation.sub.diagnosis", y = up2$gene_name[i], title = up2$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups")
#   p<- p +  geom_quasirandom()
#   plotlistup2ALLWHO[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup2ALLWHO, filename = "up2ALLWHO.pdf",
#   ncol = 3, nrow = 4
# )
# 
# rm(up2diffprotALLWHO, plotlistup2ALLWHO)
# 
# #UP3
# up3diffprotALLWHO <- pre_protein_counts_filtered2 %>% 
#   dplyr::select("patient_id", "AQP1", "FABP5") %>% 
#   dplyr::inner_join(diffprotWHO2 %>% dplyr::select(Sample_Name, methylation.sub.diagnosis), by = c("patient_id" = "Sample_Name")) %>% 
#   group_by(methylation.sub.diagnosis)
# 
# up3diffprotALLWHO$methylation.sub.diagnosis <- factor(up3diffprotALLWHO$methylation.sub.diagnosis, levels = c("A_IDH", "A_IDH_HG"))
#   
#   
# plotlistup3ALLWHO = list()
# for(i in 1:length(up3$gene_name)){
#   p <- ggboxplot(up3diffprotALLWHO, x = "methylation.sub.diagnosis", y = up3$gene_name[i], title = up3$gene_name[i],
#                  color = c("#00AFBB", "#E7B800"), ylab = "log2(Intensity)", xlab = "Groups")
#   p<- p +  geom_quasirandom()
#   plotlistup3ALLWHO[[i]]<- p
# }
# 
# ggexport (
#   plotlist = plotlistup3ALLWHO, filename = "up3ALLWHO.pdf",
#   ncol = 2, nrow = 1
# )
# 
# 
# 
# ggboxplot(up3diffprotALLWHO, x = "methylation.sub.diagnosis", y = "AQP1", add = "dotplot", color = "colors", palette =c("#00AFBB", "#E7B800"))
# 
# 



#######LIMMA NEW GRADE


# library(limma)
# 
# annotations <- readxl::read_xlsx('~/projects/glass//data//glass/Proteomics/2022-03-31_data_update/Annotation_Reduced_withControls.xlsx')
# 
# 
# diffprotWHO3 <- inner_join(diffprotWHO, WHOclass0305, by = c("Customer_ID"="Sample_Name"))
# 
# 
# y4 <- pre_protein_counts_filtered2 %>% 
#   tibble::column_to_rownames("patient_id") %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(GB_GIV_153_R1 = NULL) %>% 
#   dplyr::mutate(AC_GII_131_P = NULL) %>% 
#   dplyr::mutate(AC_GII_146_P = NULL)
# 
# 
# design4 <- model.matrix(~diffprotWHO3$Grade)
# head(design4)
# 
# 
# fit3 <- lmFit(y4, design4)
# fit3 <- eBayes(fit3)
# diffprot3 <- topTable(fit3,sort="none",n=Inf) 
# 
# save(diffprot3, file = "diffprot3.Rdata")
# 
# write.xlsx(diffprot3, file ="diffprot3.xlsx")
# 
# volcanodiffprot2 <- diffprot3 %>% 
#   rownames_to_column("gene") %>% 
#   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75,  "Significant", "Not Significant")) %>% 
#   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & gene %in% DNA_binding$name, "DNA binding", colors)) %>% 
#   dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & gene %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>% 
#   dplyr::mutate(size = ifelse(colors == "DNA binding" | colors == "Cell cycle", 1.5,1))
# 
# 
# 
# ggplot(volcanodiffprot2, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
#   geom_point(aes(col=colors, size = size)) + 
#   scale_radius(range = c(1,2))+
#   guides(size = "none")+
#   scale_color_manual(values=c("purple","royalblue", "grey", "black")) + 
#   theme_bw()+
#   geom_text_repel(data=subset(volcanodiffprot2, colors == "Cell cycle" | colors == "DNA binding"), aes(label=gene)) +
#   geom_hline(yintercept=1.3, linetype="dashed")+
#   geom_vline(xintercept = 0.75, linetype="dashed")+
#   geom_vline(xintercept = -0.75, linetype="dashed")+
#   labs(title = 'Differentially expressed genes for proteomics data', subtitle= 'Between WHO2021 grade 2&3 (n=38) vs grade 4 (n=14)')+
#   xlab('Log fold change')+
#   ylab('-log10(adjusted p-value)')
# 
# diffprot3 <- diffprot3 %>% 
#   rownames_to_column("Protein")
# 
# diffprot3signup <- diffprot3 %>% 
#   dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC >0.75)
# 
# diffprot3signup <- diffprot3signup[order(-(diffprot3signup$logFC)),]
# 
# diffprot3signdown <- diffprot3 %>% 
#   dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC < -0.75 )
# 
# diffprot3signdown <- diffprot3signdown[order((diffprot3signdown$logFC)),]
# 
# rm(fit3, y4, design4)
# 


#ENRICHMENT METHYLATION

upprot_names <- gconvert(diffprot2signup$Protein)

domain <- gconvert(corRNAprothigh$gene)

prot_up_gp_2 = prot_up_gp_2$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]



prot_up_gp = gost(list(upprot_names$name), multi_query = F, evcodes = T, sources = "GO:BP", custom_bg = domain$target)
prot_up_gp = prot_up_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
prot_up_gp <- prot_up_gp[order(prot_up_gp$p_value),]
prot_up_gp$GeneRatio = paste0(prot_up_gp$intersection_size, "/", prot_up_gp$query_size)
prot_up_gp$BgRatio = paste0(prot_up_gp$term_size, "/", prot_up_gp$effective_domain_size)
names(prot_up_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
prot_up_gp$geneID = gsub(",","/", prot_up_gp$geneID)
row.names(prot_up_gp) = prot_up_gp$ID
prot_up_gp_cluster = new("compareClusterResult", compareClusterResult = prot_up_gp)
prot_up_gp_enrich = new("enrichResult", result = prot_up_gp)
prot_up_gp_enrich@result$p.adjust<- prot_up_gp_enrich@result$p.adjust %>% 
  signif(., 2)
result_prot_up_gp <- prot_up_gp_enrich@result

barplot(prot_up_gp_enrich, showCategory = 10, fontsize = 5)


prot_up_gp_gg <-ggbarplot(result_prot_up_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")

pdf("upregprot.pdf", width = 6, height = 4)

ggpar(prot_up_gp_gg, title = "Upregulated pathways protein for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()

rm(enrichsignupprot, prot_up_gp_cluster)



#################

downprot_names <- gconvert(diffprot2signdown$Protein)

prot_down_gp = gost(list(downprot_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
prot_down_gp = prot_down_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
prot_down_gp <- prot_down_gp[order(prot_down_gp$p_value),]
prot_down_gp$GeneRatio = paste0(prot_down_gp$intersection_size, "/", prot_down_gp$query_size)
prot_down_gp$BgRatio = paste0(prot_down_gp$term_size, "/", prot_down_gp$effective_domain_size)
names(prot_down_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
prot_down_gp$geneID = gsub(",","/", prot_down_gp$geneID)
row.names(prot_down_gp) = prot_down_gp$ID
prot_down_gp_cluster = new("compareClusterResult", compareClusterResult = prot_down_gp)
prot_down_gp_enrich = new("enrichResult", result = prot_down_gp)
prot_down_gp_enrich@result$p.adjust<- prot_down_gp_enrich@result$p.adjust %>% 
  signif(., 2)
result_prot_down_gp <- prot_down_gp_enrich@result

prot_down_gp_gg <-ggbarplot(result_prot_down_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")


pdf("downregprot.pdf", width = 6.5, height = 4)

ggpar(prot_down_gp_gg, title = "Downregulated pathways protein for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()


ggplot(result_prot_down_gp[1:10,], aes (x = Count, y = Description, fill = p.adjust))+
  geom_col()+
  scale_color_gradient2(low = "black", high= "white")


save(result_prot_down_gp, file = "result_prot_down_gp.Rdata")
save(result_prot_up_gp, file = "result_prot_up_gp.Rdata")

# ###########GRADE
# 
# upprotGRADE_names <- gconvert(diffprot3signup$Protein)
# 
# prot_up_gpGRADE = gost(list(upprotGRADE_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
# prot_up_gpGRADE = prot_up_gpGRADE$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
# prot_up_gpGRADE <- prot_up_gpGRADE[order(prot_up_gpGRADE$p_value),]
# prot_up_gpGRADE$GeneRatio = paste0(prot_up_gpGRADE$intersection_size, "/", prot_up_gpGRADE$query_size)
# prot_up_gpGRADE$BgRatio = paste0(prot_up_gpGRADE$term_size, "/", prot_up_gpGRADE$effective_domain_size)
# names(prot_up_gpGRADE) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
# prot_up_gpGRADE$geneID = gsub(",","/", prot_up_gpGRADE$geneID)
# row.names(prot_up_gpGRADE) = prot_up_gpGRADE$ID
# prot_up_gpGRADE_cluster = new("compareClusterResult", compareClusterResult = prot_up_gpGRADE)
# prot_up_gpGRADE_enrich = new("enrichResult", result = prot_up_gpGRADE)
# prot_up_gpGRADE_enrich@result$p.adjust<- prot_up_gpGRADE_enrich@result$p.adjust %>% 
#   signif(., 2)
# result_prot_up_gpGRADE <- prot_up_gpGRADE_enrich@result
# 
# barplot(prot_up_gpGRADE_enrich, showCategory = 10, fontsize = 5)
# 
# 
# prot_up_gpGRADE_gg <-ggbarplot(result_prot_up_gpGRADE[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")
# 
# ggpar(prot_up_gpGRADE_gg, title = "Upregulated pathways protein between WHO2021 grade 2&3 and grade 4", font.main = c(15, "bold", "black"), xlab = "Description", ylab= "Count", 
#       font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"))
# 
# rm(enrichsignupprot, prot_up_gpGRADE_cluster)
# 
# 
# 
# #################
# 
# downprotGRADE_names <- gconvert(diffprot3signdown$Protein)
# 
# prot_down_gpGRADE = gost(list(downprotGRADE_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
# prot_down_gpGRADE = prot_down_gpGRADE$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
# prot_down_gpGRADE <- prot_down_gpGRADE[order(prot_down_gpGRADE$p_value),]
# prot_down_gpGRADE$GeneRatio = paste0(prot_down_gpGRADE$intersection_size, "/", prot_down_gpGRADE$query_size)
# prot_down_gpGRADE$BgRatio = paste0(prot_down_gpGRADE$term_size, "/", prot_down_gpGRADE$effective_domain_size)
# names(prot_down_gpGRADE) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
# prot_down_gpGRADE$geneID = gsub(",","/", prot_down_gpGRADE$geneID)
# row.names(prot_down_gpGRADE) = prot_down_gpGRADE$ID
# prot_down_gpGRADE_cluster = new("compareClusterResult", compareClusterResult = prot_down_gpGRADE)
# prot_down_gpGRADE_enrich = new("enrichResult", result = prot_down_gpGRADE)
# prot_down_gpGRADE_enrich@result$p.adjust<- prot_down_gpGRADE_enrich@result$p.adjust %>% 
#   signif(., 2)
# result_prot_down_gpGRADE <- prot_down_gpGRADE_enrich@result
# 
# prot_down_gpGRADE_gg <-ggbarplot(result_prot_down_gpGRADE[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")
# 
# 
# ggpar(prot_down_gpGRADE_gg, title = "Downregulated pathways protein between WHO2021 grade 2&3 and grade 4", font.main = c(15, "bold", "black"), xlab = "Description", ylab= "Count", 
#       font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"))
# 
# 
# 
# ggplot(result_prot_down_gpGRADE[1:10,], aes (x = Count, y = Description, fill = p.adjust))+
#   geom_col()+
#   scale_color_gradient2(low = "black", high= "white")
# 
# #https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
# 
# 
# 
#-----------------Sample_Type
table(diffprotWHO$methylation.sub.diagnosis)
diffprotWHO$methylation.sub.diagnosis

y3 <- pre_protein_counts_filtered2 %>%
  tibble::column_to_rownames("patient_id") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(GB_GIV_153_R1 = NULL) %>%
  dplyr::mutate(AC_GII_131_P = NULL) %>%
  dplyr::mutate(AC_GII_146_P = NULL) %>%
  dplyr::mutate(AC_GII_136_P = NULL) %>%
  dplyr::mutate(AC_GII_141_R1 = NULL)


diffprotWHO2 <- diffprotWHO %>%
  dplyr::filter(methylation.sub.diagnosis %in% c("O_IDH", "CONTR_HEMI") == F) %>%
  dplyr::filter(Sample_Name %in% c("AC_GII_131_P", "AC_GII_146_P", "GB_GIV_153_R1") ==F)

diffprotWHO2_paired <- diffprotWHO2 %>%
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Customer_ID)) %>%
  dplyr::mutate(batch = ifelse(Sample_Number %in% c("105","111","113","115","117","122","128","137","146","148","150","156","157","174"), "001", Sample_Number)) %>% 
  dplyr::mutate(batch = paste0("s",batch)) %>% 
  dplyr::mutate(batch = factor(batch)) %>% 
  dplyr::mutate(methylation.sub.diagnosis  = factor(methylation.sub.diagnosis , levels=c('A_IDH','A_IDH_HG'))) 
  
stopifnot(length(diffprotWHO2_paired$batch %>% levels) == length(diffprotWHO2_paired$Sample_Number[duplicated(diffprotWHO2_paired$Sample_Number)])+1)



design3_paired <- model.matrix(~diffprotWHO2_paired$batch +diffprotWHO2_paired$methylation.sub.diagnosis)
head(design3_paired)


fit_paired <- lmFit(y3, design3_paired)
fit_paired <- eBayes(fit_paired)
diffprot2_paired <- topTable(fit_paired,sort="none",n=Inf)

diffprot2_paired <- diffprot2_paired %>% 
  dplyr::select(diffprotWHO2_paired.methylation.sub.diagnosisA_IDH_HG, AveExpr, `F`, P.Value, adj.P.Val) %>% 
  rename("logFC" = diffprotWHO2_paired.methylation.sub.diagnosisA_IDH_HG) %>% 
  rownames_to_column("Protein")


diffprot2_paired_signup <- diffprot2_paired %>%
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC >0.75)

diffprot2_paired_signup <- diffprot2_paired_signup[order(-(diffprot2_paired_signup$logFC)),]

diffprot2_paired_signdown <- diffprot2_paired %>%
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC < -0.75 )

diffprot2_paired_signdown <- diffprot2_paired_paired_signdown[order((diffprot2_paired_paired_signdown$logFC)),]

#patprot

# table(diffprotWHO$methylation.sub.diagnosis)
# diffprotWHO$methylation.sub.diagnosis
# 
# y3 <- pre_protein_counts_filtered2 %>%
#   tibble::column_to_rownames("patient_id") %>%
#   t() %>%
#   as.data.frame() %>%
#   dplyr::mutate(GB_GIV_153_R1 = NULL) %>%
#   dplyr::mutate(AC_GII_131_P = NULL) %>%
#   dplyr::mutate(AC_GII_146_P = NULL) %>%
#   dplyr::mutate(AC_GII_136_P = NULL) %>%
#   dplyr::mutate(AC_GII_141_R1 = NULL)
# 
# 
# 
# table(metadatapatprot$methylation.sub.diagnosis)
# metadatapatprot$methylation.sub.diagnosis
# 
# 
# design3_paired <- model.matrix(~diffprotWHO2_paired$Sample_Number + diffprotWHO2_paired$methylation.sub.diagnosis)
# head(design3_paired)
# 
# 
# fit_paired <- lmFit(y3, design3_paired)
# fit_paired <- eBayes(fit_paired)
# diffprot2_paired <- topTable(fit_paired,sort="none",n=Inf)

diffprotWHO$Sample_Number


diffprotWHO$Sample_Number[!duplicated(diffprotWHO$Sample_Number)]

duplicates_prot <- diffprotWHO$Sample_Number[duplicated(diffprotWHO$Sample_Number)]

sample_unique_prot <- diffprotWHO %>% 
  dplyr::filter(Sample_Number %!in% duplicates_prot)
length(sample_unique_prot$Customer_ID)



## LGHG

y3 <- pre_protein_counts_filtered2 %>% 
  tibble::column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(GB_GIV_153_R1 = NULL) %>% 
  dplyr::mutate(AC_GII_131_P = NULL) %>% 
  dplyr::mutate(AC_GII_146_P = NULL) %>% 
  dplyr::mutate(AC_GII_136_P = NULL) %>% 
  dplyr::mutate(AC_GII_141_R1 = NULL)


diffprotWHOLGHG <- diffprotWHO %>% 
  dplyr::filter(methylation.sub.diagnosis %in% c("O_IDH", "CONTR_HEMI") == F) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_131_P", "AC_GII_146_P", "GB_GIV_153_R1") ==F) %>% 
  dplyr::select(Customer_ID, methylation.sub.diagnosis, Sample_Name, Sample_Number) %>% 
  dplyr::filter(!(Customer_ID == "104_R2" | Customer_ID == "164_R1" | Customer_ID == "102_R2" | Customer_ID == "103_R2" |
                    Customer_ID == "119_R2" | Customer_ID == "124_R2" | Customer_ID == "143_R2" | Customer_ID == "126_R2"
                  | Customer_ID == "108_R2" | Customer_ID == "130_R2" | Customer_ID == "139_R2"))




y3LGHG <- y3 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample_Name") %>% 
  dplyr::filter(Sample_Name %in% diffprotWHOLGHG$Sample_Name) %>% 
  column_to_rownames("Sample_Name") %>% 
  t() %>% 
  as.data.frame()

design3LGHG <- model.matrix(~diffprotWHOLGHG$methylation.sub.diagnosis)
head(design3LGHG)


fitLGHG <- lmFit(y3LGHG, design3LGHG)
fitLGHG <- eBayes(fitLGHG)
diffprot2LGHG <- topTable(fitLGHG,sort="none",n=Inf) 

diffprot2LGHG <- diffprot2LGHG %>% 
  rownames_to_column("Protein")

diffprot2LGHG <- diffprot2LGHG[order(diffprot2LGHG$adj.P.Val),]



diffprot2signupLGHG <- diffprot2LGHG %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC >0.75)

diffprot2signupLGHG <- diffprot2signupLGHG[order(-(diffprot2signupLGHG$logFC)),]

diffprot2signdownLGHG <- diffprot2LGHG %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC < -0.75 )


diffprotWHOLGHG$Sample_Number[duplicated(diffprotWHOLGHG$Sample_Number)]

ion_transport <- read_csv("ion_transport.csv")


volcanodiffprotLGHG <- diffprot2LGHG %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75,  "Significant", "Not Significant")) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & logFC > 0.75 & Protein %in% DNA_replication$name, "DNA replication", colors)) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & Protein %in% cycling.cell.markers$gene_name, "Cell cycle", colors)) %>%
  dplyr::mutate(colors =ifelse(adj.P.Val <0.05 & abs(logFC) > 0.75 & Protein %in% ion_transport$name, "Ion transport", colors)) %>%
  dplyr::mutate(size = ifelse(colors == "DNA replication" | colors == "Cell cycle" | colors == "Ion transport", 1.5,1))



pdf("Volcano_prot.pdf", width = 10, height = 7)
ggplot(volcanodiffprotLGHG, aes(logFC, -log10(adj.P.Val), color = colors, label = Protein)) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(data = subset(volcanodiffprotLGHG, colors == "Not Significant"), size = 1)+
  geom_point(data = subset(volcanodiffprotLGHG, colors == "Significant"), size = 1)+
  geom_point(data = subset(volcanodiffprotLGHG, colors == "DNA replication"), size = 2)+
  geom_point(data = subset(volcanodiffprotLGHG, colors == "Cell cycle"), size = 2)+
  geom_point(data = subset(volcanodiffprotLGHG, colors == "Ion transport"), size = 2)+
  scale_color_manual(values = c("gray87", "gray40", "royalblue", "purple", "red"),
                     name = "colors",
                     breaks=c("Not Significant", "Significant", "DNA replication", "Cell cycle", "Ion transport"),
                     labels = c("Not Significant", "Significant", "DNA replication", "Cell cycle", "Ion transport"))+
  theme_bw()+
  geom_text_repel(data=subset(volcanodiffprotLGHG, colors == "Cell cycle" | colors == "DNA replication" | colors == "Ion transport"), aes(label=Protein), color = "black") +
  geom_hline(yintercept=1.3, linetype="dashed")+
  geom_vline(xintercept = 0.75, linetype="dashed")+
  geom_vline(xintercept = -0.75, linetype="dashed")+
  labs(title = 'Differentially expressed genes for proteomics data', subtitle= 'Between methylation classifier LG (n=29) versus HG (n=10)')+
  xlab('Log fold change')+
  ylab('-log10(adjusted p-value)')+
  xlim(-2.3,2.3)


dev.off()


#ENRICHMENT METHYLATION

upprot_names <- gconvert(diffprot2signupLGHG$Protein)


prot_up_gp = gost(list(upprot_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
prot_up_gp = prot_up_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
prot_up_gp <- prot_up_gp[order(prot_up_gp$p_value),]
prot_up_gp$GeneRatio = paste0(prot_up_gp$intersection_size, "/", prot_up_gp$query_size)
prot_up_gp$BgRatio = paste0(prot_up_gp$term_size, "/", prot_up_gp$effective_domain_size)
names(prot_up_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
prot_up_gp$geneID = gsub(",","/", prot_up_gp$geneID)
row.names(prot_up_gp) = prot_up_gp$ID
prot_up_gp_cluster = new("compareClusterResult", compareClusterResult = prot_up_gp)
prot_up_gp_enrich = new("enrichResult", result = prot_up_gp)
prot_up_gp_enrich@result$p.adjust<- prot_up_gp_enrich@result$p.adjust %>% 
  signif(., 2)
result_prot_up_gp <- prot_up_gp_enrich@result

barplot(prot_up_gp_enrich, showCategory = 10, fontsize = 5)


prot_up_gp_gg <-ggbarplot(result_prot_up_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkblue")

pdf("upregprot.pdf", width = 7, height = 4)

ggpar(prot_up_gp_gg, title = "Upregulated pathways protein for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()

rm(enrichsignupprot, prot_up_gp_cluster)



#################

downprot_names <- gconvert(diffprot2signdownLGHG$Protein)

prot_down_gp = gost(list(downprot_names$name), multi_query = F, evcodes = T, sources = "GO:BP")
prot_down_gp = prot_down_gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
prot_down_gp <- prot_down_gp[order(prot_down_gp$p_value),]
prot_down_gp$GeneRatio = paste0(prot_down_gp$intersection_size, "/", prot_down_gp$query_size)
prot_down_gp$BgRatio = paste0(prot_down_gp$term_size, "/", prot_down_gp$effective_domain_size)
names(prot_down_gp) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
prot_down_gp$geneID = gsub(",","/", prot_down_gp$geneID)
row.names(prot_down_gp) = prot_down_gp$ID
prot_down_gp_cluster = new("compareClusterResult", compareClusterResult = prot_down_gp)
prot_down_gp_enrich = new("enrichResult", result = prot_down_gp)
prot_down_gp_enrich@result$p.adjust<- prot_down_gp_enrich@result$p.adjust %>% 
  signif(., 2)
result_prot_down_gp <- prot_down_gp_enrich@result

prot_down_gp_gg <-ggbarplot(result_prot_down_gp[1:10,], "Description", "Count", orientation = "horiz", fill = "darkred")


pdf("downregprot.pdf", width = 6.5, height = 4)

ggpar(prot_down_gp_gg, title = "Downregulated pathways protein for high-grade tumors", font.main = c(9, "bold", "black"), xlab = "Description", ylab= "Count", 
      font.x = c(8, "plain", "black"), font.y = c(8, "plain", "black"), font.ytickslab = c(7, "plain", "black"), font.xtickslab =  7)

dev.off()




##ONLY 1 SAMPLE
diffprotWHOLGHG$Sample_Number[duplicated(diffprotWHOLGHG$Sample_Number)]

y3LGHG <- y3 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample_Name") %>% 
  dplyr::filter(Sample_Name %in% diffprotWHOLGHG$Sample_Name) %>% 
  column_to_rownames("Sample_Name") %>% 
  t() %>% 
  as.data.frame()

design3LGHG <- model.matrix(~diffprotWHOLGHG$methylation.sub.diagnosis)
head(design3LGHG)


fitLGHG <- lmFit(y3LGHG, design3LGHG)
fitLGHG <- eBayes(fitLGHG)
diffprot2LGHG <- topTable(fitLGHG,sort="none",n=Inf) 

diffprot2LGHG <- diffprot2LGHG %>% 
  rownames_to_column("Protein")

diffprot2LGHG <- diffprot2LGHG[order(diffprot2LGHG$adj.P.Val),]

diffprot2LGHG$Protein[diffprot2LGHG$Protein == "TARSL2"] <-"TARS3"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "ASNA1"] <-"GET3"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "EPRS"] <-"EPRS1"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "ADSS"] <-"ADSS2"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "CRAD"] <-"CRACD"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "HIST1H2BA"] <-"H2BC1"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "KIF1BP"] <-"KIFBP"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "FAM129B"] <-"NIBAN2"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "FAM49A"] <-"CYRIA"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "RYDEN"] <-"SHFL"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "FAM49B"] <-"CYRIB"
diffprot2LGHG$Protein[diffprot2LGHG$Protein == "ADPRHL2"] <-"ADPRS"



diffprot2signupLGHG <- diffprot2LGHG %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC >0.75)

diffprot2signupLGHG <- diffprot2signupLGHG[order(-(diffprot2signupLGHG$logFC)),]

diffprot2signdownLGHG <- diffprot2LGHG %>% 
  dplyr::filter(adj.P.Val <as.numeric(0.05) & logFC < -0.75 )
