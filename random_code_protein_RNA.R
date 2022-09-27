#lijnplots RNA en protein ----
first20genesRNA <- total_RNA[order(total_RNA$gene_name),] %>% 
  dplyr::slice_head(n=10) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("gene_name")

first20genesprotein <- protein_counts_filtered[order(protein_counts_filtered$Protein),] %>% 
  dplyr::slice_head(n=10) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("Protein")

par(mfrow = c(2,1))

matplot(t(first20genesRNA), type = "l")
legend("topleft", legend = rownames(first20genesRNA), col=1:10, lty=1:10, cex=0.6)

matplot(t(first20genesprotein), type = "l")
legend("topleft", legend = rownames(first20genesprotein), col=1:10, lty=1:10, cex=0.6)
dev.off()

##plots per gen correlatie----
A1Bgprot <- t(first20genesprotein[1,]) %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id")

A1Bgprot<-A1Bgprot[order(A1Bgprot$patient_id),] %>% 
  as.data.frame() %>% 
  column_to_rownames("patient_id")

A1Bgprot <- as.data.frame(column_to_rownames("patient_id"))

A1BgRNA = t(first20genesRNA[1,]) %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id")

A1BgRNA<-A1BgRNA[order(A1BgRNA$patient_id),] 

totalA1Bg <-left_join(A1BgRNA, A1Bgprot, by ="patient_id")

plot(totalA1Bg$A1BG.x, totalA1Bg$A1BG.y)
corA1Bg <- cor(as.numeric(totalA1Bg$A1BG.x), as.numeric(totalA1Bg$A1BG.y))   


#from the 200 most abundant genes, which ones are shared between RNA and protein
mostabundantprot <- highprotexpr2 %>% 
  slice_head(n=200) %>% 
  tibble::rownames_to_column("protein")

mostabundantRNA <- highRNAexpr2 %>% 
  slice_head(n=200) %>% 
  tibble::rownames_to_column("RNA")

sharedmostabundant <- inner_join(mostabundantprot, mostabundantRNA, by =c("protein"="RNA"))

rm(mostabundantprot, mostabundantRNA)

#100 proteins with highest protein correlation 
corr_check <- function(Dataset, threshold){
  matriz_cor <- cor(Dataset)
  matriz_cor
  
  for (i in 1:nrow(matriz_cor)){
    correlations <-  which((abs(matriz_cor[i,i:ncol(matriz_cor)]) > threshold) & (matriz_cor[i,i:ncol(matriz_cor)] != 1))
    
    if(length(correlations)> 0){
      print(correlations)
      
    }
  }
}

corr_check(t(highprotexpr2), 0.9)

library(lares)



#PCA proteomics----
pca_prot <- mostabundant_prot %>% 
  #tibble: column_to_rownames("Protein") #moet erbij als tmpprot wordt gebruikt
  as.data.frame() %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('Sample_Name')

pca_total <- dplyr::left_join(metadata_prot_RNA,pca_prot, by=c("Sample_Name.x" = "Sample_Name"))

ggplot(pca_total, aes(x=PC1, y=PC2, color = Sample_Sex)) +
  geom_point() +
  theme_classic()


#intensities before normalization
protein <- tmpprot %>% 
  tibble::column_to_rownames("Protein")

boxplot(protein)


#intensities after normalization
normalized_protein <- log2(protein)
boxplot(normalized_protein)

#summarized experiment----
#t = make_se(d %>%  dplyr::rename(name = X) %>%  dplyr::mutate(ID = name), which(colnames(d) != 'X'), metadata_prot_RNA)

proteins_unique = data.frame(tmpprot[,1])

exprs <- tmpprot %>% 
  tibble::column_to_rownames('Protein')
colData <- data.frame(label = metadata_protRNA$Customer_ID, condition = metadata_prot_RNA$resec2, replicate = replicate(length(metadata_prot_RNA$resec2),0))
rowRanges <- rowRanges(tmpprot[,-1])
SE = SummarizedExperiment(assay = list("exprs"=exprs),colData=colData, metadata = metadata_prot_RNA)

data_diff <- test_diff(SE)


columns <- 1:length(colData$label);columns
SE2 = make_se(exprs, columns,colData)

tmpprot_unique <- tmpprot %>% 
  dplyr::mutate(ids = c(1:length(tmpprot[,1])))

proteins_unique <- tmpprot [,1]