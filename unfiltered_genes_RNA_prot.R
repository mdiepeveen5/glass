#needsGS_ID_RNA_proteinand protein_counts_filtered

readcounts_unfiltered <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(readcounts_unfiltered) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(readcounts_unfiltered)))

readcounts_unfiltered2<- readcounts_unfiltered %>%  
  dplyr::mutate(Chr=NULL) %>% 
  dplyr::mutate(Start=NULL) %>% 
  dplyr::mutate(Length=NULL) %>% 
  dplyr::mutate(End=NULL) %>% 
  dplyr::mutate(Strand=NULL) %>% 
  dplyr::filter(Geneid %in% c("ENSG00000280987.4" , "ENSG00000187522.16","ENSG00000285053.1", "ENSG00000285441.1", "ENSG00000033050.9", "ENSG00000169093.16_PAR_Y", "ENSG00000169100.14_PAR_Y") == F) %>% 
  dplyr::mutate(gene_id = NULL) %>%
  tibble::column_to_rownames("Geneid") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('GS_ID') %>% 
  dplyr::mutate("GS_ID" = gsub(".","-", GS_ID, fixed = T)) %>% 
  dplyr::right_join(GS_ID_RNA_protein) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene_name") %>% 
  dplyr::inner_join(gene.annot2, by= c("gene_name" = "gene_id")) %>% 
  dplyr::mutate(gene_name=NULL) %>% 
  dplyr::inner_join(dplyr::select(protein_counts_filtered_allprotein, Protein), by= c("gene_name.y" = "Protein")) %>% 
  dplyr::relocate("gene_name.y") %>% 
  dplyr::rename("gene_name" = "gene_name.y")


interprotRNAunfiltered <- dplyr::inner_join(readcounts_unfiltered2, dplyr::select(protein_counts_filtered_allprotein, Protein), by= c("gene_name" = "Protein")) %>% 
  dplyr::select(gene_name)

readcounts_unfiltered2<-readcounts_unfiltered2 %>% 
  column_to_rownames("gene_name")


metadata2<-
  data.frame(GS_ID = colnames(readcounts_unfiltered2)) %>% 
  dplyr::left_join(
    read.csv("~/projects/glass/data/glass/Clinical data/Cleaned/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") ,
    by=c('GS_ID'='GS_ID')
  ) %>% 
  dplyr::mutate(Exclude=NULL) %>%   
  dplyr::mutate(Surgery_ID=NULL) %>% 
  dplyr::mutate(Sample_Sex = ifelse(is.na(Sample_Sex),"unknown",Sample_Sex))

dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts_unfiltered2,
                                       colData = metadata2,
                                       design= ~Sample_Sex)
dds2 <-DESeq(dds2)
vsd2 <- DESeq2::vst(dds2, blind=T)
readcounts_unfiltered_vst <- assay(vsd2) %>% as.data.frame 

rm(vsd2, dds2, readcounts_unfiltered, metadata2)


proteinordered_unfiltered <- dplyr::left_join(interprotRNAunfiltered, protein_counts_filtered_allprotein, by = c("gene_name" = "Protein")) %>% 
  dplyr::select(c('gene_name', sortpatientid))


proteinordered_unfiltered <- proteinordered_unfiltered[order(proteinordered_unfiltered$gene_name),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL) 


readcounts_unfiltered_vst_ordered <- readcounts_unfiltered_vst %>% 
  rownames_to_column("gene_name") %>% 
  dplyr::select(c('gene_name', sortpatientid))

readcounts_unfiltered_vst_ordered<- readcounts_unfiltered_vst_ordered[order(readcounts_unfiltered_vst_ordered$gene_name),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL)  

rm(readcounts_unfiltered, readcounts_unfiltered2)

corproteinordered_unfiltered <- proteinordered_unfiltered %>% 
  column_to_rownames('gene_name') %>% 
  t() %>% 
  as.data.frame()

rm(proteinordered_unfiltered)

corRNAordered_unfiltered <- readcounts_unfiltered_vst_ordered %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame()

rm(readcounts_unfiltered_vst_ordered)

corRNAprotOLDSPEAR_unfiltered=data.frame()
for(i in colnames(corRNAordered_unfiltered)[1:length(corRNAordered_unfiltered)]) {
  print(i)
  a <- corRNAordered_unfiltered %>% dplyr::select(all_of(i))
  b <- corproteinordered_unfiltered %>% dplyr::select(all_of(i))
  
  
  corRNAprotOLDSPEAR_unfiltered <- rbind(corRNAprotOLDSPEAR_unfiltered, data.frame(gene=i, 
                                                             correlation = as.numeric(cor(a,b, method = "spearman"))
                                                    
  ))
}
rm(i,a,b)


corRNAprotOLDSPEAR_unfiltered <- corRNAprotOLDSPEAR_unfiltered[order(-corRNAprotOLDSPEAR_unfiltered$correlation),]

highRNAexpr_unfiltered <- t(corRNAordered_unfiltered) %>% 
  as.data.frame() %>% 
  dplyr::mutate(sumrow = rowSums(.))
highRNAexpr_unfiltered <- highRNAexpr_unfiltered[order(-highRNAexpr_unfiltered$sumrow),] 

highRNAexprplot <- highRNAexpr_unfiltered %>% 
  dplyr::select(sumrow) %>% 
  tibble::rownames_to_column("gene_name")

highRNAcor_unfiltered <- left_join(highRNAexprplot, corRNAprotOLDSPEAR_unfiltered, by = c("gene_name" = "gene")) 

highRNAcor_unfiltered$group[highRNAcor_unfiltered$gene_name %in% c(interprotRNA$gene_name)] <- "Unfiltered"
highRNAcor_unfiltered$group[!(highRNAcor_unfiltered$gene_name %in% c(interprotRNA$gene_name))] <- "Filtered"


ggplot(data = highRNAcor_unfiltered, aes (x = sumrow, y = correlation, color = group))+
  geom_point()+
  geom_smooth(data = subset(highRNAcor_unfiltered, group != 'Filtered')) +
  labs ()+
  xlab ("Sumrow of RNA expression per gene")+
  ylab ("Correlation Protein x RNA")

cor.test(log( 1+highRNAcor$sumrow), highRNAcor$correlation)

rm(highRNAexprplot, highRNAcor)


rm(interprotRNAunfiltered, highRNAexpr_unfiltered, highRNAcor_unfiltered, tmp, unfiltered_RNA, unfiltered_protein, readcounts_unfiltered_vst)


