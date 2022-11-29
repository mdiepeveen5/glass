load("annotations.Rdata", "pre_protein_counts_filered2")



#purity is bepaald op DNA-seq 

methylation.purities <- read.table("data/glass/Methylation/Analysis/RFpurity/purities_RFpurity.txt") %>% 
  dplyr::mutate(methylation.uid = paste0(Slide, "_", Array)) %>% 
  dplyr::mutate(fn = NULL) %>% 
  dplyr::mutate(Basename = NULL) %>% 
  dplyr::mutate(Slide = NULL) %>% 
  dplyr::mutate(Array = NULL) %>% 
  dplyr::rename(methylation.purity.absolute = absolute) %>% 
  dplyr::rename(methylation.purity.estimate = estimate)

stopifnot(duplicated(methylation.purities$methylation.uid) == F) # no duplicates may exist



tmp.methylation.metdata <- read.csv("data/glass/Methylation/Metadata/Datasheet4.csv") %>% 
  dplyr::mutate(methylation.uid = paste0(Slide, "_", Array)) %>% 
  dplyr::mutate(
    Basename = NULL,
    Array = NULL,
    Slide = NULL,
    Surgery_ID = NULL, 
    GLASS_ID = NULL, 
    Sample_Plate = NULL,       
    Sample_ID = NULL, 
    Sample_Resection = NULL, 
    Sample_Type = NULL, 
    Recurrent_Type = NULL, 
    Sample_Sex = NULL  )

stopifnot(duplicated(tmp.methylation.metdata$methylation.uid) == F) # no duplicates may exist


# I ran all 323 samples I found on the server with RFpurity, only a subset matches the actual glass samples (tmp.2)
stopifnot(tmp.methylation.metdata$methylation.uid %in% methylation.purities$methylation.uid)

methylation.purities <-methylation.purities %>% 
  dplyr::left_join(tmp.methylation.metdata, by=c('methylation.uid'='methylation.uid'))

rm(tmp.methylation.metdata)


## append to metadata ----


VAFpurities <- read.csv("purities_2.csv") %>% 
  dplyr::select("Sample_Name", "dna.wes.VAF_IDH")



metadataALL <- metadataALL %>% 
  dplyr::left_join(methylation.purities, by=c('Sample_Name'='Sample_Name'), keep=F,suffix = c("", "")) %>% 
  dplyr::left_join(VAFpurities, by = c("Sample_Name"))# force overwrite



metapurity<- inner_join(diffprotWHO, metadataALL, by = c("Customer_ID"="Sample_Name", "methylation.sub.diagnosis")) %>% 
  dplyr::select("Sample_Name", "methylation.purity.estimate", "methylation.sub.diagnosis") %>% 
  dplyr::filter(Sample_Name %in% c("GB_GIV_147_P", "GB_GIV_147_R2", "AC_GII_131_P","GB_GIV_103_R1", "AC_GII_126_R1", "AC_GII_146_P", "AC_GII_130_R1", "GB_GIV_153_R1") == F) %>% 
  left_join(pre_protein_counts_filtered2 %>% dplyr::select("patient_id","PLP1", "PLLP", "MAG","MOG"), by = c("Sample_Name"="patient_id"))




ggplot(metapurity, aes(y = methylation.purity.estimate))+
  geom_point(aes(x = PLP1), color = "red1")+
  geom_point(aes(x = PLLP), color = "dodgerblue3")+
  geom_point(aes(x = MAG), color = "darkolivegreen")+
  geom_point(aes(x = MOG), color = "purple4")+
  theme_bw()+ 
  xlab("Protein")+
  ylab ("Tumor Purity")
  


cor.test(metapurity$methylation.purity.estimate, metapurity$PLP1)
#-0.13 pval 0.37
cor.test(metapurity$methylation.purity.estimate, metapurity$PLLP)
#-0.09 pval 0.56
cor.test(metapurity$methylation.purity.estimate, metapurity$MAG)
#-0.11 pval 0.48
cor.test(metapurity$methylation.purity.estimate, metapurity$MOG)
#-0.14 pval 0.36


t.test (metapurity %>% dplyr::filter(methylation.sub.diagnosis == "A_IDH") %>% dplyr::select(methylation.purity.estimate), metapurity %>% dplyr::filter(methylation.sub.diagnosis == "A_IDH_HG") %>% dplyr::select(methylation.purity.estimate), var.equal = T)

VAFpurity <- inner_join(diffprotWHO, metadataALL, by = c("Customer_ID"="Sample_Name", "methylation.sub.diagnosis")) %>% 
  dplyr::select("Sample_Name", "dna.wes.VAF_IDH", "methylation.sub.diagnosis") %>% 
  rename("VAF" = "dna.wes.VAF_IDH") %>% 
  dplyr::filter(Sample_Name %in% c("AC_GII_136_P") == F) %>% 
  left_join(pre_protein_counts_filtered2 %>% dplyr::select("patient_id","PLP1", "PLLP", "MAG","MOG"), by = c("Sample_Name"="patient_id")) %>% 
  na.omit()


ggplot(VAFpurity, aes(y = VAF))+
  geom_point(aes(x = PLP1), color = "red1")+
  geom_point(aes(x = PLLP), color = "dodgerblue3")+
  geom_point(aes(x = MAG), color = "darkolivegreen")+
  geom_point(aes(x = MOG), color = "purple4")+
  theme_bw()+ 
  xlab("Protein")+
  ylab ("Tumor Purity")

rm(VAFpurity)

cor.test(VAFpurity$VAF, VAFpurity$PLP1)
cor.test(VAFpurity$VAF, VAFpurity$PLLP)
cor.test(VAFpurity$VAF, VAFpurity$MAG)
cor.test(VAFpurity$VAF, VAFpurity$MOG)

t.test (VAFpurity %>% dplyr::filter(methylation.sub.diagnosis == "A_IDH") %>% dplyr::select(VAF), VAFpurity %>% dplyr::filter(methylation.sub.diagnosis == "A_IDH_HG") %>% dplyr::select(VAF), var.equal = T)
