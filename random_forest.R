library(randomForestSRC)
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

tmp.1 <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Surgery data_GLASS RNAseq.csv')
tmp.1 <- rbind(
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S1")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.)))
  ,
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S2")) %>% 
    `colnames<-`(gsub('_S[1-4]$','',colnames(.))),
  tmp.1 %>%
    dplyr::select(`GLASS_ID` | ends_with("_S3")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.))),
  tmp.1 %>% dplyr::select(`GLASS_ID` | ends_with("_S4")) %>%
    `colnames<-`(gsub('_S[1-4]$','',colnames(.)))
) %>% 
  tidyr::drop_na(Sample_Name) %>% 
  dplyr::arrange(Sample_Name)
#dplyr::select(Sample_Name, Date_Surgery, GLASS_ID)


tmp.2 <- read.csv('data/glass/Clinical data/Cleaned/metadata_2022/Survival data_GLASS RNAseq__ALL.csv') %>% 
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Diagnosis , format = "%Y-%m-%d")) %>% 
  dplyr::mutate(Date_of_Diagnosis = as.Date(Date_of_Death , format = "%Y-%m-%d")) %>% 
  dplyr::select(GLASS_ID, Date_of_Death, Date_Last_Followup) %>% 
  dplyr::mutate(Date_Last_Event = ifelse(is.na(Date_of_Death), Date_Last_Followup , Date_of_Death)) %>% 
  dplyr::mutate(Date_Last_Event.status = ifelse(is.na(Date_of_Death), 0 , 1) ) %>% 
  dplyr::mutate(Date_of_Death = NULL,  Date_Last_Followup = NULL)


stopifnot(tmp.1$GLASS_ID %in% tmp.2$GLASS_ID)


tmp <- tmp.1 %>% dplyr::left_join(tmp.2, by=c('GLASS_ID'='GLASS_ID')) %>% 
  dplyr::mutate(time.resection.until.last.event = difftime(Date_Last_Event,Date_Surgery, units = 'days')) %>% 
  dplyr::rename(status.resection.until.last.event = Date_Last_Event.status) %>% 
  dplyr::select(Sample_Name, time.resection.until.last.event, status.resection.until.last.event)


stopifnot(tmp$time.resection.until.last.event > 0)
rm(tmp.1, tmp.2)



metadataDE2 <- metadataDE %>%
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name'),suffix = c("", ""))


#metadata.glass.per.resection %>%  dplyr::filter(is.na(status.resection.until.last.event))


rm(tmp)

#works if load Change_names_and_order
#needs metadataALL, sortpatientid

#152 genes 
MADprotein <- corproteinordered %>%
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)

# plot(sort(MADprotein$mad, decreasing = T))


MADRNA <- corRNAordered %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)


# plot(sort(MADRNA$mad, decreasing = T))

#ONLY TO MAKE 'OVERLAP' 
MADRNA2 <- MADRNA %>% 
  dplyr::filter(mad>1)%>% 
  dplyr::mutate(mad =NULL) %>% 
  rownames_to_column("gene")


MADprotein2 <- MADprotein %>% 
  dplyr::filter(mad>1)%>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("protein")

overlap <- union(MADRNA2$gene, MADprotein2$protein)

rm(MADRNA2, MADprotein2)


MADRNAsurv <- MADRNA %>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% overlap) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID")




survival_RNArfc <- annotations %>% dplyr::select(GS_ID, Deceased, OS_FirstSurgery_months, Survival_FirstRecurrenceSurgery_months, Condition) %>%
  dplyr::mutate(NewSurvival = ifelse(Condition == "P", OS_FirstSurgery_months, Survival_FirstRecurrenceSurgery_months)) %>% 
  dplyr::filter(Condition != "R2") %>% 
  dplyr::inner_join(MADRNAsurv) %>% 
  dplyr::mutate(GS_ID = NULL) %>% 
  dplyr::mutate(OS_FirstSurgery_months = NULL) %>% 
  dplyr::mutate(Survival_FirstRecurrenceSurgery_months = NULL) %>% 
  dplyr::mutate(Condition = NULL)
  


  
survival_RNA_plot <- rfsrc(Surv(NewSurvival,Deceased) ~ ., survival_RNArfc , importance = TRUE)
plot.rfsrc(survival_RNA_plot)

boruta.survival <- Boruta(survival_RNA_firstrecur~., data = survival_RNA_firstrecur)


plot.survival(survival_RNA_firstrecurplot)

#PROTEIN
MADprotsurv <- MADprotein %>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% overlap) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID")


survival_protrfc <- annotations %>% dplyr::select(GS_ID, Deceased, Survival_FirstRecurrenceSurgery_months, Condition) %>% 
  dplyr::inner_join(MADprotsurv) %>% 
  dplyr::mutate(GS_ID = NULL)


survival_prot_firstrecur <- survival_protrfc %>% 
  dplyr::filter(Condition == "R") %>% 
  dplyr::mutate(Condition = NULL)


survival_prot_firstrecurplot <- rfsrc(Surv(Survival_FirstRecurrenceSurgery_months,Deceased) ~ ., survival_prot_firstrecur , importance = TRUE)
plot.rfsrc(survival_prot_firstrecurplot)

plot.survival(survival_prot_firstrecurplot)


#MAD after first recurrence
GS_ID_firstrecur <- annotations %>% 
  dplyr::filter(Condition == "R")

MADproteinfirstrecur <- corproteinordered %>%
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% c(GS_ID_firstrecur$GS_ID)) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)

plot(sort(MADproteinfirstrecur$mad, decreasing = T))


MADproteinfirstrecur2 <- MADproteinfirstrecur %>% 
  dplyr::slice_head(n= 10) %>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("protein")

#RNA
MADRNAfirstrecur <- corRNAordered %>%
  rownames_to_column("GS_ID") %>% 
  dplyr::filter(GS_ID %in% c(GS_ID_firstrecur$GS_ID)) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad)

plot(sort(MADRNAfirstrecur$mad, decreasing = T))


MADRNAfirstrecur2 <- MADRNAfirstrecur %>% 
  dplyr::slice_head(n= 10) %>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("gene")


overlap2 <- union(MADRNAfirstrecur2$gene, MADproteinfirstrecur2$protein)

MADRNAsurv2 <- MADRNAfirstrecur %>% 
  dplyr::mutate(mad = NULL) %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% overlap) %>% 
  column_to_rownames("gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID")


survival_RNArfc2 <- annotations %>% dplyr::select(GS_ID, Deceased, Survival_FirstRecurrenceSurgery_months) %>% 
  dplyr::inner_join(MADRNAsurv2) %>% 
  dplyr::mutate(GS_ID = NULL)



survival_RNA_firstrecurplot2 <- rfsrc(Surv(Survival_FirstRecurrenceSurgery_months,Deceased) ~ ., survival_RNArfc2 , importance = TRUE)
plot.rfsrc(survival_RNA_firstrecurplot2)

plot.survival(survival_RNA_firstrecurplot2)



#MAD gesoorteerd plotten en kijken waar hij afloopt
#1e en dan 2e en dan allebei


#How many patients low grade high grade after first recurrence in last dataset:
# patients <- metadataDE %>% dplyr::filter(GS_ID %in% sortpatientid) %>% 
#   left_join(annotations) %>% 
#   dplyr::filter(Condition == "R")



data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 3, importance = T)

plt <- data.frame(importance = pbc.obj$importance) %>% 
  dplyr::arrange(abs(importance)) %>% 
  dplyr::top_n(20) %>% 
  dplyr::arrange(-importance) %>% 
  tibble::rownames_to_column('gene_uid') %>% 
  dplyr::mutate(baseline = 0) %>% 
  dplyr::mutate(y = rank(importance)) %>% 
  tidyr::pivot_longer(cols = -c(gene_uid, y)) %>% 
  dplyr::mutate(name=NULL) %>% 
  dplyr::rename(RFSRC.importance.R2.exon = value)


plot.rfsrc(pbc.obj)
# plot.variable.rfsrc(pbc.obj)
plot.survival(pbc.obj)



#tijd van surgery tot dood gebruiken 


MADRNAall <- readcount.vstDE %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::filter(mad > 1) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate("gene_id" = NULL) %>% 
  select(gene_name, everything())%>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("GS_ID")



survival_RNAallrfc <- metadataDE2 %>% dplyr::select(GS_ID, status.resection.until.last.event, time.resection.until.last.event, resection) %>% 
  dplyr::filter(resection == "S2") %>% 
  dplyr::mutate(resection = NULL) %>% 
  dplyr::inner_join(MADRNAall) %>% 
  dplyr::mutate(GS_ID = NULL) %>% 
  na.omit() %>% 
  transform(time.resection.until.last.event = as.numeric(time.resection.until.last.event))
            
  



survivalRNAall_plot <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~ ., survival_RNAallrfc , importance = TRUE)
plot.rfsrc(survivalRNAall_plot)

plot.survival(survivalRNAall_plot)

#proteomics
MADprotall <- pre_protein_counts_filtered2%>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::filter(mad > 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id") %>% 
  left_join(diffprotWHO %>% dplyr::select(Sample_Name, Customer_ID), by = c("patient_id" = "Sample_Name")) %>% 
  dplyr::mutate(patient_id = NULL) %>%
  select(Customer_ID, everything())%>% 
  na.omit() 

survival_protallrfc <- annotations %>% 
  dplyr::select(Sample_Name, Deceased, Survival_FirstRecurrenceSurgery_months, Condition) %>% 
  dplyr::filter(Condition == "P") %>% 
  dplyr::mutate(Condition = NULL) %>% 
  dplyr::inner_join(MADprotall, by = c("Sample_Name" = "Customer_ID")) %>% 
  dplyr::mutate(Sample_Name = NULL) %>% 
  na.omit()



survival_protallrfc_plot <- rfsrc(Surv(Survival_FirstRecurrenceSurgery_months,Deceased) ~ ., survival_protallrfc)
plot.rfsrc(survival_protallrfc_plot)

plot.survival(survival_protallrfc_plot)

print(survival_protallrfc_plot)
############resec2
survival_protallrfcR2 <- annotations %>% 
  dplyr::select(Sample_Name, Deceased, Survival_FirstRecurrenceSurgery_months, Condition) %>% 
  dplyr::filter(Condition == "R") %>% 
  dplyr::mutate(Condition = NULL) %>% 
  dplyr::inner_join(MADprotall, by = c("Sample_Name" = "Customer_ID")) %>% 
  dplyr::mutate(Sample_Name = NULL) %>% 
  na.omit()



survival_protallrfc_plotR2 <- rfsrc(Surv(Survival_FirstRecurrenceSurgery_months,Deceased) ~ ., survival_protallrfcR2)
plot.rfsrc(survival_protallrfc_plotR2)

plot.survival(survival_protallrfc_plotR2)

print(survival_protallrfc_plotR2)

#######################
MADprotall2 <- MADprotall %>% 
  inner_join(tmp, by = c("Customer_ID" = "Sample_Name")) %>% 
  dplyr::mutate(Customer_ID = NULL) %>% 
  transform(time.resection.until.last.event = as.numeric(time.resection.until.last.event)) %>% 
  select(time.resection.until.last.event, everything())%>% 
  select(status.resection.until.last.event, everything()) %>% 
  dplyr::mutate(status.resection.until.last.event = 1)
  
#stopifnot(length(MADprotall2$batch %>% levels) == length(MADprotall2$Sample_Number[duplicated(MADprotall2$Sample_Number)])+1)


#importance uit 

survival_protallrfc_plot <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., MADprotall2 , ntree = 25)
plot.rfsrc(survival_protallrfc_plot, sorted = T, verbose = T, cex = 0.5)

print(survival_protallrfc_plot)
plot(gg_vimp(survival_protallrfc_plot), nvar =20)

plot.survival(survival_protallrfc_plot)

summary(survival_protallrfc_plot)

print(survival_protallrfc_plot)

##MAD OMLAAG GEEFT SLECHTERE PREDICTIE

survival_protallrfc_plot.pred <- predict(survival_protallrfc_plot, MADprotall2, outcome = "test")


data (pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status)~. , pbc, importance = T)
plot.survival(pbc.obj)
plot.rfsrc(pbc.obj)
plot(get.brier.survival(pbc.obj, cens.model = "rfsrc")$brier.score, type = "s", col =4)

print(pbc.obj)


MADRNA2 <- readcount.vst %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::filter(mad > 1) %>%
  dplyr::mutate(mad= NULL) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::select(gene_name, everything()) %>% 
  dplyr::mutate(gene_id = NULL) %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id") %>% 
  left_join(metadata %>% dplyr::select(Sample_Name, GS_ID), by = c("patient_id" = "GS_ID")) %>% 
  dplyr::mutate(patient_id = NULL) %>%
  select(Sample_Name, everything())%>% 
  na.omit() %>% 
  inner_join(tmp, by = c("Sample_Name" = "Sample_Name")) %>% 
  dplyr::mutate(Sample_Name =NULL) %>% 
  transform(time.resection.until.last.event = as.numeric(time.resection.until.last.event)) %>% 
  select(time.resection.until.last.event, everything())%>% 
  select(status.resection.until.last.event, everything()) 

survival_RNAallrfc_plot <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., MADRNA2 , importance = T)
plot.rfsrc(survival_RNAallrfc_plot, sorted = T, verbose = T, cex = 0.5)

plot(gg_vimp(survival_RNAallrfc_plot), nvar =40)


print(survival_RNAallrfc_plot)


#####################################################

MADRNAsurv <- readcount.vst %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::filter(mad > 1) %>%
  dplyr::mutate(mad= NULL) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::select(gene_name, everything()) %>% 
  dplyr::mutate(gene_id = NULL) %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id") %>% 
  left_join(metadata %>% dplyr::select(Sample_Name, GS_ID, resection), by = c("patient_id" = "GS_ID")) %>% 
  dplyr::mutate(patient_id = NULL) %>%
  select(Sample_Name, everything())%>% 
  na.omit() %>% 
  inner_join(tmp, by = c("Sample_Name" = "Sample_Name")) %>% 
  transform(time.resection.until.last.event = as.numeric(time.resection.until.last.event)) %>% 
  select(time.resection.until.last.event, everything())%>% 
  select(status.resection.until.last.event, everything()) %>% 
  select(resection, everything())

MADRNAsurvR1 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S1") %>% 
  dplyr::mutate(resection = NULL) %>% 
  dplyr::mutate(Sample_Name = NULL)

survival_RNAallrfc_plotR1 <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., MADRNAsurvR1 , importance = "permute")
plot.rfsrc(survival_RNAallrfc_plotR1, sorted = T, verbose = T, cex = 0.5)


survival_RNAallrfc_plotR1$var.used

VIMPsR1<- as.data.frame(survival_RNAallrfc_plotR1$importance) %>% 
  rename("Importance" = "survival_RNAallrfc_plotR1$importance")


plot(gg_vimp(survival_RNAallrfc_plotR1), nvar =20, main = "Hallo")

print(survival_RNAallrfc_plotR1)



MADRNAsurvR2 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S2") %>% 
  dplyr::mutate(resection = NULL) %>% 
  dplyr::mutate(Sample_Name = NULL)

survival_RNAallrfc_plotR2 <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., MADRNAsurvR2 , importance = "permute")
plot.rfsrc(survival_RNAallrfc_plotR2, sorted = T, verbose = T, cex = 0.5)

plot(gg_vimp(survival_RNAallrfc_plotR2), nvar =40)

print(survival_RNAallrfc_plotR2)



VIMPsR2<- as.data.frame(survival_RNAallrfc_plotR2$importance) %>% 
  rename("Importance" = "survival_RNAallrfc_plotR2$importance") %>% 
  rownames_to_column("Gene") %>% 
  dplyr::mutate(colors = ifelse(Gene %in% cycling.cell.markers$gene_name, "Cell Cycle", "Other"))

VIMPsR2 <- VIMPsR2[order(-VIMPsR2$Importance),] 

rownames(VIMPsR2) <- NULL


pdf("VIMP.pdf", height = 4, width = 8)
ggbarplot(VIMPsR2[1:10,], "Gene", "Importance", color = "colors", palette = c("purple", "red"))
dev.off()
#PCNA los coxph


get.cindex(time = MADRNAsurvR2$time.resection.until.last.event, censoring = MADRNAsurvR2$status.resection.until.last.event, predicted = survival_RNAallrfc_plotR2$predicted.oob)

###RNA with less patients


patients_R2 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S2") %>% 
  dplyr::mutate(resection = NULL) %>%
  dplyr::select(Sample_Name) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*", "", Sample_Name)) %>% 
  dplyr::mutate(Sample_Name = NULL)

patients_R1 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S1") %>% 
  dplyr::mutate(resection = NULL) %>%
  dplyr::select(Sample_Name) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*", "", Sample_Name)) %>% 
  dplyr::mutate(Sample_Name = NULL)



lesspatientsR2 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S2") %>% 
  dplyr::mutate(resection = NULL) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*", "", Sample_Name)) %>% 
  dplyr::mutate(Sample_Name = NULL) %>% 
  dplyr::filter(Sample_Number %in% patients_R1$Sample_Number) %>% 
  dplyr::mutate(Sample_Number = NULL)


bothdatasets <- patients_R2 %>% 
  anti_join(lesspatientsR2, by = c("Sample_Number"))



survival_RNAallrfc_plotR2_lesspatients <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., lesspatientsR2 , importance = "permute")
plot.rfsrc(survival_RNAallrfc_plotR2_lesspatients, sorted = T, verbose = T, cex = 0.5)

plot(gg_vimp(survival_RNAallrfc_plotR2_lesspatients), nvar =40)

print(survival_RNAallrfc_plotR2_lesspatients)


#R1

lesspatientsR1 <- MADRNAsurv %>% 
  dplyr::filter(resection == "S1") %>% 
  dplyr::mutate(resection = NULL) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*", "", Sample_Name)) %>% 
  dplyr::mutate(Sample_Name = NULL) %>% 
  dplyr::filter(Sample_Number %in% patients_R2$Sample_Number) %>% 
  dplyr::mutate(Sample_Number = NULL)
  
survival_RNAallrfc_plotR1_lesspatients <- rfsrc(Surv(time.resection.until.last.event,status.resection.until.last.event) ~., lesspatientsR1 , importance = TRUE)
plot.rfsrc(survival_RNAallrfc_plotR1_lesspatients, sorted = T, verbose = T, cex = 0.5)

plot(gg_vimp(survival_RNAallrfc_plotR1_lesspatients), nvar =40)

print(survival_RNAallrfc_plotR1_lesspatients)
