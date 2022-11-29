tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))


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
                                      design=~1)
dds <-DESeq(dds)
vsd <- DESeq2::vst(dds)

pca_unsupervised <- assay(vsd) %>% as.data.frame %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=1000) %>% 
  dplyr::mutate(mad = NULL) %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')

unsupervised_sorted <-dplyr::right_join(metadata, pca_unsupervised) %>% 
  dplyr::mutate(ID = str_sub(GLASS_ID, -3))


##PCA plot----

pdf("PCA_RNA.pdf", height = 6,width = 10)
ggplot(unsupervised_sorted%>%  dplyr::arrange(Sample_Type, Sample_Name), aes(x=PC1, y=PC2, color = ID, group = ID)) +
  geom_point() +
  geom_path(arrow=arrow(), type="closed") +
  theme_classic()

dev.off()
differencesPCA <- as.data.frame(unsupervised_sorted %>% dplyr::select(GLASS_ID, PC1, PC2)) %>% 
  dplyr::filter(GLASS_ID %in% patient_sums) %>% 
  dplyr::mutate(PC1_rand = sample(PC1)) %>% 
  dplyr::mutate(PC2_rand = sample(PC2))

sumsPCA = data.frame()
for (i in 1:length(differencesPCA$GLASS_ID)){
  
  if (differencesPCA$GLASS_ID[i] == differencesPCA$GLASS_ID[i+1]){
    sumsPCA <- rbind(sumsPCA, data.frame(GLASS_ID = differencesPCA$GLASS_ID[i], differencePCA1 = abs(differencesPCA$PC1[i] - differencesPCA$PC1[i+1]),
                                         differencePCA2 = abs(differencesPCA$PC2[i] - differencesPCA$PC2[i+1]),
                                         differencePCA1_rand = abs(differencesPCA$PC1[i] - differencesPCA$PC1_rand[i+1]),
                                         differencePCA2_rand = abs(differencesPCA$PC1[i] - differencesPCA$PC2_rand[i+1])
                                         ))
  }
} 

patient_sums <- sumsPCA$GLASS_ID


mean(sumsPCA$differencePCA1)
mean(sumsPCA$differencePCA1_rand)
mean(sumsPCA$differencePCA2)
mean(sumsPCA$differencePCA2_rand)

wilcox.test(sumsPCA$differencePCA1, sumsPCA$differencePCA1_rand, var.equal = T)
t.test(sumsPCA$differencePCA2, sumsPCA$differencePCA2_rand, var.equal = T)


df = data.frame( a = c(sumsPCA$differencePCA2, sumsPCA$differencePCA2_rand)) %>% 
  dplyr::mutate(b = c(
    rep("A",length(sumsPCA$differencePCA1))
    ,
    rep("B",length(sumsPCA$differencePCA1_rand))
  ))

ggplot(df, aes(x=b, y=a)) + geom_jitter()

#sumsPCA %>% dplyr::filter(abs(differencePCA1) < 0.1*(max(differencesPCA$PC1) - min(differencesPCA$PC1)))

#sumsPCA %>% dplyr::filter(abs(differencePCA2) < 0.1*(max(differencesPCA$PC2) - min(differencesPCA$PC2)))


rm(dds, vsd, res, selected.features.unsupervised, pca_supervised, unsupervised_sorted)




library(limma)



y5 <- pre_protein_counts_filtered2 %>%
  tibble::column_to_rownames("patient_id") %>%
  t() %>%
  as.data.frame() %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=1000) %>% 
  dplyr::mutate(mad = NULL) 



pca_unsupervised_prot <- y5 %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('Sample_Name')

unsupervised_sorted_prot <-dplyr::inner_join(diffprotWHO %>% dplyr::select(Customer_ID, Sample_Name), pca_unsupervised_prot) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Customer_ID)) %>% 
  dplyr::mutate(Sample_Type = gsub(".*_", "", Customer_ID))

##PCA plot----
pdf("PCA_prot.pdf", height = 6, width = 9)
ggplot(unsupervised_sorted_prot%>%  dplyr::arrange(Sample_Type, Customer_ID), aes(x=PC1, y=PC2, color = Sample_Number, group = Sample_Number)) +
  geom_point() +
  geom_path(arrow=arrow(), type="closed") +
  theme_classic()+
  labs(color = "ID")


dev.off()

differencesPCA_prot <- as.data.frame(unsupervised_sorted_prot %>% dplyr::select(Sample_Number, PC1, PC2)) %>% 
  dplyr::filter(Sample_Number %in% patient_sums_prot) %>% 
  dplyr::mutate(PC1_rand = sample(PC1)) %>% 
  dplyr::mutate(PC2_rand = sample(PC2))

sumsPCA_prot = data.frame()
for (i in 1:length(differencesPCA_prot$Sample_Number)){
  
  if (differencesPCA_prot$Sample_Number[i] == differencesPCA_prot$Sample_Number[i+1]){
    sumsPCA_prot <- rbind(sumsPCA_prot, data.frame(Sample_Number = differencesPCA_prot$Sample_Number[i], differencePCA1 = abs(differencesPCA_prot$PC1[i] - differencesPCA_prot$PC1[i+1]),
                                         differencePCA2 = abs(differencesPCA_prot$PC2[i] - differencesPCA_prot$PC2[i+1]),
                                         differencePCA1_rand = abs(differencesPCA$PC1[i] - differencesPCA$PC1_rand[i+1]),
                                         differencePCA2_rand = abs(differencesPCA$PC1[i] - differencesPCA$PC2_rand[i+1])))
  }
} 

patient_sums_prot <- sumsPCA_prot$Sample_Number

mean(sumsPCA_prot$differencePCA1)
mean(sumsPCA_prot$differencePCA1_rand)
mean(sumsPCA_prot$differencePCA2)
mean(sumsPCA_prot$differencePCA2_rand)

t.test(sumsPCA_prot$differencePCA1, sumsPCA_prot$differencePCA1_rand, var.equal = T)
t.test(sumsPCA_prot$differencePCA2, sumsPCA_prot$differencePCA2_rand, var.equal = T)


sumsPCA_prot %>% dplyr::filter(abs(differencePCA1) < 0.1*(max(differencesPCA_prot$PC1) - min(differencesPCA_prot$PC1)))


df = data.frame( a = c(sumsPCA_prot$differencePCA1,sumsPCA_prot$differencePCA1_rand)) %>% 
  dplyr::mutate(b = c(
    rep("A",length(sumsPCA_prot$differencePCA1))
    ,
    rep("B",length(sumsPCA_prot$differencePCA1_rand))
  ))

ggplot(df, aes(x=b, y=a)) + geom_jitter()

df = data.frame( a = c(sumsPCA_prot$differencePCA2, sumsPCA_prot$differencePCA2_rand)) %>% 
  dplyr::mutate(b = c(
    rep("A",length(sumsPCA_prot$differencePCA1))
    ,
    rep("B",length(sumsPCA_prot$differencePCA1_rand))
  ))

ggplot(df, aes(x=b, y=a)) + geom_jitter()


sumsPCA_prot %>% dplyr::filter(abs(differencePCA2) < 0.1*(max(differencesPCA_prot$PC2) - min(differencesPCA_prot$PC2)))




##LESS META


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
readcount.vst_lastres <- assay(vsd_lastres) %>% as.data.frame 


pca_unsupervised_lastres <- assay(vsd_lastres) %>% as.data.frame %>% 
  dplyr::mutate(mad = apply(as.matrix(.), 1, stats::mad)) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=1000) %>% 
  dplyr::mutate(mad = NULL) %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')

unsupervised_sorted_lastres <-dplyr::right_join(metadata_onlylastresec, pca_unsupervised_lastres) %>% 
  dplyr::mutate(ID = str_sub(GLASS_ID, -3))

pdf("PCA_RNA.pdf", height = 6,width = 10)
ggplot(unsupervised_sorted_lastres%>%  dplyr::arrange(Sample_Type, Sample_Name), aes(x=PC1, y=PC2, color = ID, group = ID)) +
  geom_point() +
  geom_path(arrow=arrow(), type="closed") +
  theme_classic()
dev.off()

differencesPCA_lastres <- as.data.frame(unsupervised_sorted_lastres %>% dplyr::select(GLASS_ID, PC1, PC2)) %>% 
  dplyr::mutate(PC1_rand = sample(PC1)) %>% 
  dplyr::mutate(PC2_rand = sample(PC2))

sumsPCA_lastres = data.frame()
for (i in 1:length(differencesPCA_lastres$GLASS_ID)){
  
  if (differencesPCA_lastres$GLASS_ID[i] == differencesPCA_lastres$GLASS_ID[i+1]){
    sumsPCA_lastres <- rbind(sumsPCA_lastres, data.frame(GLASS_ID = differencesPCA_lastres$GLASS_ID[i], differencePCA1 = abs(differencesPCA_lastres$PC1[i] - differencesPCA_lastres$PC1[i+1]),
                                         differencePCA2 = abs(differencesPCA_lastres$PC2[i] - differencesPCA_lastres$PC2[i+1]),
                                         differencePCA1_rand = abs(differencesPCA_lastres$PC1[i] - differencesPCA_lastres$PC1_rand[i+1]),
                                         differencePCA2_rand = abs(differencesPCA_lastres$PC1[i] - differencesPCA_lastres$PC2_rand[i+1])
    ))
  }
} 


wilcox.test(sumsPCA_lastres$differencePCA1, sumsPCA_lastres$differencePCA1_rand, var.equal = T)
t.test(sumsPCA_lastres$differencePCA2, sumsPCA_lastres$differencePCA2_rand, var.equal = T)

mean(sumsPCA_lastres$differencePCA2_rand)


##PROT ON METH


unsupervised_sorted_prot_meth <-dplyr::inner_join(diffprotWHO2 %>% dplyr::select(Customer_ID, Sample_Name, methylation.sub.diagnosis), pca_unsupervised_prot) %>% 
  dplyr::mutate(Sample_Number = gsub("\\_.*","",Customer_ID)) %>% 
  dplyr::mutate(Sample_Type = gsub(".*_", "", Customer_ID))

##PCA plot----
pdf("PCA_prot.pdf", height = 6, width = 9)
ggplot(unsupervised_sorted_prot_meth%>%  dplyr::arrange(Sample_Type, Customer_ID), aes(x=PC1, y=PC2, color = methylation.sub.diagnosis, group = Sample_Number)) +
  geom_point() +
  geom_path(arrow=arrow(), type="closed") +
  theme_classic()+
  labs(color = "ID")

