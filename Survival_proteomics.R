library(dplyr)
library(tibble)
library(corrplot)
library(ggfortify)
library(survminer)

save(pre_protein_counts_filtered2, file = "pre_protein_counts_filtered2.Rdata")
load("diffprotWHO2.Rdata", "pre_protein_counts_filered2")

save(precorsignup, file = "precorsignup.Rdata")

precorsignup <- pre_protein_counts_filtered2 %>% 
  right_join(diffprotWHO2 %>% dplyr::select("Sample_Name"), by = c("patient_id" = "Sample_Name"))

corsignupprotein <- precorsignup %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  inner_join(diffprot2signup %>% dplyr::select(Protein), by = c("gene" = "Protein")) %>% 
  column_to_rownames("gene")

rm(precorsignup)


corrplot(cor(as.data.frame(t(corsignupprotein)), use = "pairwise.complete.obs", method = "spearman"), order = "hclust", tl.cex = 0.4)

?corrplot

###########################



############
annotations <- readxl::read_xlsx('~/projects/glass//data//glass/Proteomics/2022-03-31_data_update/Annotation_Reduced_withControls.xlsx')


annotations <- annotations %>% 
  dplyr::mutate("Grade" = NULL) %>% 
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

#save(annotations, file = "annotations.Rdata")

annotations_filtered <- inner_join(diffprotWHO2, annotations, by = c("Customer_ID"="Sample_Name")) %>% 
  left_join(WHOclass0305 %>% dplyr::select(Sample_Name, Grade), by = c("Customer_ID"="Sample_Name"))


annotations_filteredR2 <- annotations_filtered %>% 
  dplyr::filter(Condition == "R") 

#ALLGENES

allgenes <- corsignupprotein %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')

PCAallgenes <- allgenes %>% 
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sorted_allgenes <- inner_join(annotations_filteredR2, PCAallgenes)

PCAdefallgenes <- prcomp(allgenes, scale = T)

autoplot(PCAdefallgenes, data = PCA_sorted_allgenes, colour = "methylation.sub.diagnosis")

screeplot(PCAdefallgenes)  

Survallgenes <- PCA_sorted_allgenes %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 5, "High", "Low"))

plot(sort(PCA_sorted_allgenes$PC1))


survivalallgenes <- survfit(Surv(Survallgenes$Survival_FirstRecurrenceSurgery_months, Survallgenes$Deceased)~Type, data = Survallgenes)

ggsurvplot(survivalallgenes, pval = TRUE)

#pval = 0.0074




WHOsurvivalallgenes <- annotations_filteredR2 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, Grade) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(Grade == 4, "High", "Low"))


survivalWHO <- survfit(Surv(WHOsurvivalallgenes$Survival_FirstRecurrenceSurgery_months, WHOsurvivalallgenes$Deceased)~Type, data = WHOsurvivalallgenes)


ggsurvplot(survivalWHO, pval = TRUE)

##########CLUST1

clust1 <- corsignupprotein %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% c("RAVER1", "TMPO", "NFIX", "RCC2", "FIP1L1", "RPF2", "DDX47")) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')


PCAclust1 <- clust1 %>%  
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sortedclust1 <- inner_join(annotations_filteredR2, PCAclust1)


PCAdefclust1 <- prcomp(clust1, scale = T)

autoplot(PCAdefclust1, data = PCA_sortedclust1, colour = "methylation.sub.diagnosis")

screeplot(PCAdefclust1)  

Survclust1 <- PCA_sortedclust1 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 0, "High", "Low"))
  


survivalclust1 <- survfit(Surv(Survclust1$Survival_FirstRecurrenceSurgery_months, Survclust1$Deceased)~Type, data = Survclust1)

ggsurvplot(survivalclust1, pval = TRUE)



#pval is 0.035

#########CLUSTER2

clust2 <- corsignupprotein %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% c("COL18A1", "COL6A2", "COL6A3", "TGFBI", "BGN")) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')


PCAclust2 <- clust2 %>%  
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sortedclust2 <- inner_join(annotations_filteredR2, PCAclust2)


PCAdefclust2 <- prcomp(clust2, scale = T)


autoplot(PCAdefclust2, data = PCA_sortedclust2, colour = "methylation.sub.diagnosis")

screeplot(PCAdefclust2)  

Survclust2 <- PCA_sortedclust2 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 0, "High", "Low"))

plot(sort(PCA_sortedclust2$PC1))


survivalclust2 <- survfit(Surv(Survclust2$Survival_FirstRecurrenceSurgery_months, Survclust2$Deceased)~Type, data = Survclust2)

ggsurvplot(survivalclust2, pval = TRUE)
#pval is 0.018


##CLUSTER3

clust3 <- corsignupprotein %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% c("RFC4", "MCM3", "MCM7", "NCBP1", "H3-3A", "RCC1")) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')


PCAclust3 <- clust3 %>%  
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sortedclust3 <- inner_join(annotations_filteredR2, PCAclust3)


PCAdefclust3 <- prcomp(clust3, scale = T)


autoplot(PCAdefclust3, data = PCA_sortedclust3, colour = "methylation.sub.diagnosis")

screeplot(PCAdefclust3)  

Survclust3 <- PCA_sortedclust3 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 0, "High", "Low"))



survivalclust3 <- survfit(Surv(Survclust3$Survival_FirstRecurrenceSurgery_months, Survclust3$Deceased)~Type, data = Survclust3)

ggsurvplot(survivalclust3, pval = TRUE)
#0.047

#cluster4

clust4 <- corsignupprotein %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% c("PCNA")) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')


PCAclust4 <- clust4 %>%  
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sortedclust4 <- inner_join(annotations_filteredR2, PCAclust4)


PCAdefclust4 <- prcomp(clust4, scale = T)


autoplot(PCAdefclust4, data = PCA_sortedclust4, colour = "methylation.sub.diagnosis")

screeplot(PCAdefclust4)  

Survclust4 <- PCA_sortedclust4 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 0, "High", "Low"))



survivalclust4 <- survfit(Surv(Survclust4$Survival_FirstRecurrenceSurgery_months, Survclust4$Deceased)~Type, data = Survclust4)

ggsurvplot(survivalclust4, pval = TRUE)
#0.038
rm(clust4, PCAclust4, PCAdefclust4, PCA_sortedclust4)

########################


corsigndownprotein <- precorsignup %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  inner_join(diffprot2signdown %>% dplyr::select(Protein), by = c("gene" = "Protein")) %>% 
  column_to_rownames("gene")



clust5 <- corsigndownprotein %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% c("MBP")) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')


PCAclust5 <- clust5 %>%  
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sortedclust5 <- inner_join(annotations_filteredR2, PCAclust5)


PCAdefclust5 <- prcomp(clust5, scale = T)


autoplot(PCAdefclust5, data = PCA_sortedclust5, colour = "methylation.sub.diagnosis")

screeplot(PCAdefclust5)  

Survclust5 <- PCA_sortedclust5 %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(status = ifelse(Deceased == 1, 2,1)) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 0, "High", "Low"))



survivalclust5 <- survfit(Surv(Survclust5$Survival_FirstRecurrenceSurgery_months, Survclust5$status)~Type, data = Survclust5)

ggsurvplot(survivalclust5, pval = TRUE)


#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/



precorsignup2 <- pre_protein_counts_filtered2 %>% 
  right_join(diffprotWHO3 %>% dplyr::select("Sample_Name"), by = c("patient_id" = "Sample_Name"))

corsignupprotein2 <- precorsignup2 %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  inner_join(diffprot2signup %>% dplyr::select(Protein), by = c("gene" = "Protein")) %>% 
  column_to_rownames("gene")

rm(precorsignup2)

annotations_filtered3 <- inner_join(diffprotWHO3, annotations, by = c("Customer_ID"="Sample_Name"))


annotations_filteredR2_3 <- annotations_filtered3 %>% 
  dplyr::filter(Condition.x == "R") 



allgenes2 <- corsignupprotein2 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2_3$Sample_Name) %>% 
  column_to_rownames('Sample_Name')

PCAallgenes2 <- allgenes2 %>% 
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sorted_allgenes2 <- inner_join(annotations_filteredR2_3, PCAallgenes2)


PCAdefallgenes2 <- prcomp(allgenes2, scale = T)

autoplot(PCAdefallgenes2, data = PCA_sorted_allgenes2, colour = "newgrade")

screeplot(PCAdefallgenes2)  

plot(sort(PCA_sorted_allgenes2$PC1))

abline(h = 0)

Survallgenes2 <- PCA_sorted_allgenes2 %>% 
  dplyr::select(Sample_Name, Customer_ID, Deceased.x, Survival_FirstRecurrenceSurgery_months.x, PC1) %>%
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  transform(Survival_FirstRecurrenceSurgery_months.x = as.numeric(Survival_FirstRecurrenceSurgery_months.x),
            Deceased.x = as.numeric(Deceased.x)) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 4, "High", "Low"))



survivalallgenes2 <- survfit(Surv(Survallgenes2$Survival_FirstRecurrenceSurgery_months.x, Survallgenes2$Deceased.x)~Type, data = Survallgenes2)

ggsurvplot(survivalallgenes2, pval = TRUE)

#pval = 0.0015


precorsignupGRADE <- pre_protein_counts_filtered2 %>% 
  right_join(diffprotWHO3 %>% dplyr::select("Sample_Name"), by = c("patient_id" = "Sample_Name"))



corsignupproteinGRADE <- precorsignupGRADE %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  inner_join(diffprot3signup %>% dplyr::select(Protein), by = c("gene" = "Protein")) %>% 
  column_to_rownames("gene")

rm(precorsignupGRADE)

allgenesGRADE <- corsignupproteinGRADE %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample_Name') %>% 
  dplyr::filter(Sample_Name %in% annotations_filteredR2$Sample_Name) %>% 
  column_to_rownames('Sample_Name')

PCAallgenesGRADE <- allgenesGRADE %>% 
  prcomp() %>% 
  purrr::pluck('x') %>%
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')  

PCA_sorted_allgenesGRADE <- inner_join(annotations_filteredR2, PCAallgenesGRADE)

PCAdefallgenesGRADE <- prcomp(allgenesGRADE, scale = T)


screeplot(PCAdefallgenesGRADE)  

SurvallgenesGRADE <- PCA_sorted_allgenesGRADE %>% 
  dplyr::select(Sample_Name, Sample_Condition, Customer_ID, Deceased, Survival_FirstRecurrenceSurgery_months, PC1) %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_130_R1")== F) %>% 
  dplyr::mutate(Type = ifelse(PC1 > 5, "High", "Low"))

plot(sort(PCA_sorted_allgenes$PC1))


survivalallgenesGRADE <- survfit(Surv(SurvallgenesGRADE$Survival_FirstRecurrenceSurgery_months, SurvallgenesGRADE$Deceased)~Type, data = SurvallgenesGRADE)

ggsurvplot(survivalallgenesGRADE, pval = TRUE)



###########################

coxsurvival <- rfsrc(Surv(Survallgenes$Survival_FirstRecurrenceSurgery_months, Survallgenes$Deceased)~Type, data = Survallgenes)


survivalallgenes <- survfit(Surv(Survallgenes$Survival_FirstRecurrenceSurgery_months, Survallgenes$Deceased)~Type, data = Survallgenes)


#############

save(pre_protein_counts_filtered2, file = "pre_protein_counts_filtered2.Rds")

PCNA <- pre_protein_counts_filtered2 %>% 
  dplyr::select(patient_id,PCNA) %>% 
  left_join(annotations_filtered %>% dplyr::select(Sample_Name, OS_FirstSurgery_months,Survival_FirstRecurrenceSurgery_months, Deceased, Condition), by = c("patient_id"="Sample_Name")) %>% 
  drop_na

PCNAR1 <- PCNA %>% 
  dplyr::filter(Condition == "P") %>% 
  dplyr::mutate(Type = ifelse(PCNA >6, "High", "Low"))





survivalPCAR1 <- survfit(Surv(PCNAR1$OS_FirstSurgery_months, PCNAR1$Deceased)~Type, data = PCNAR1)

ggsurvplot(survivalPCAR1, pval = TRUE)




PCNAR2 <- PCNA %>% 
  dplyr::filter(Condition == "R") %>% 
  dplyr::mutate(Type = ifelse(PCNA >6, "High", "Low"))


survivalPCAR2 <- survfit(Surv(PCNAR2$Survival_FirstRecurrenceSurgery_months, PCNAR2$Deceased)~Type, data = PCNAR2)

ggsurvplot(survivalPCAR2, pval = TRUE)



survdiff(Surv(PCNAR2$Survival_FirstRecurrenceSurgery_months, PCNAR2$Deceased)~Type, data = PCNAR2)

str(survdiff(Surv(PCNAR2$Survival_FirstRecurrenceSurgery_months, PCNAR2$Deceased)~Type, data = PCNAR2))


#coxph met hazard ratios