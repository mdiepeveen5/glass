load("NanoString_protein_data.RDS")
load("NanoString_protein_metadata.RDS")

load("NanoString_RNA_data.RDS")
load("NanoString_RNA_metadata.RDS")



NanoString_protein_data_col <- as.data.frame(NanoString_protein_data) %>% 
  tibble::rownames_to_column("Protein")
  
NanoString_RNA_data_col<- as.data.frame(NanoString_RNA_data) %>% 
  tibble::rownames_to_column("Gene")



NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Fibronectin"] <-"FN1"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD45"] <-"PTPRC"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Bcl-2"] <-"BCL2"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "SMA"] <-"SMN1"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD66b"] <-"CEACAM8"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "NY-ESO-1"] <-"CTAG1B"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "PR"] <-"PGR"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Her2"] <-"ERBB2"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Tim-3"] <-"HAVCR2"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "B7-H3"] <-"CD276"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "VISTA"] <-"VISR"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "GITR"] <-"TNFRSF18"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD25"] <-"IL2RA"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD127"] <-"IL7R"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD25"] <-"IL2RA"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Ki-67"] <-"MKI67"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD11c"] <-"ITGAX"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "PD-L1"] <-"CD274"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "Beta-2-microglobulin"] <-"B2M"
NanoString_protein_data_col$Protein[NanoString_protein_data_col$Protein == "CD8"] <-"CD8A"

spatial_prot_rna <- inner_join(NanoString_protein_data_col, NanoString_RNA_data_col, by = c("Protein" = "Gene"))


notinrnaprot <- anti_join(NanoString_protein_data_col, NanoString_RNA_data_col, by = c("Protein" = "Gene"))



NanoString_protein_data_col_filtered <- NanoString_protein_data_col %>% 
  dplyr::filter(Protein %in% spatial_prot_rna$Protein)

  
  
NanoString_protein_data_col_filtered <- NanoString_protein_data_col_filtered[order(NanoString_protein_data_col_filtered$Protein),]
  

NanoString_protein_data_col_filtered <- NanoString_protein_data_col_filtered %>% 
  dplyr::mutate(E205_010 = NULL) %>% 
  dplyr::mutate(F206_001 = NULL) %>% 
  dplyr::mutate(F206_002 = NULL) %>% 
  as.data.frame(row.names = 1:nrow(.)) %>% 
  tibble::column_to_rownames("Protein") %>% 
  dplyr::select(sort(current_vars())) %>% 
  t() %>% 
  as.data.frame()



NanoString_RNA_data_col_filtered<- NanoString_RNA_data_col %>% 
  dplyr::filter(Gene %in% spatial_prot_rna$Protein)%>% 
  dplyr::mutate(C203_012 = NULL) %>% 
  dplyr::mutate(B202_011 = NULL) %>% 
  dplyr::mutate(B202_012 = NULL) 

NanoString_RNA_data_col_filtered <- NanoString_RNA_data_col_filtered[order(NanoString_RNA_data_col_filtered$Gene),]


NanoString_RNA_data_col_filtered <- NanoString_RNA_data_col_filtered %>% 
  tibble::column_to_rownames("Gene") %>% 
  dplyr::select(sort(current_vars())) %>% 
  t() %>% 
  as.data.frame()


rm(NanoString_protein_data, NanoString_protein_data_col, NanoString_RNA_data, NanoString_RNA_data_col)



SpatialCOR=data.frame()
for(i in colnames(NanoString_protein_data_col_filtered)[1:length(NanoString_protein_data_col_filtered)]) {
  print(i)
  a <- NanoString_protein_data_col_filtered %>% dplyr::select(all_of(i))
  b <- NanoString_RNA_data_col_filtered %>% dplyr::select(all_of(i))
  
  
  print(paste0("Gene ",i,": ",cor(a,b)))
  
  SpatialCOR <- rbind(SpatialCOR, data.frame(gene=i, 
                                            correlation = as.numeric(cor(a,b, method = "spearman")),
                                            correlation.rnd = as.numeric(cor(permute::shuffle(a),b, method = "spearman"))
  
  ))
}


SpatialCOR <- SpatialCOR[order(-SpatialCOR$correlation),]


SpatialCORcompare <- SpatialCOR %>%
  dplyr::mutate(order1 = rank(-correlation, ties.method = "first")) %>% 
  dplyr::mutate(order2 = rank(-correlation.rnd, ties.method = "first")) 

Spatcorrna1 <- SpatialCORcompare %>% 
  dplyr::select(gene, correlation,order1) %>% 
  dplyr::mutate(Status = "Correlation")

Spatcorrna2 <- SpatialCORcompare %>% 
  dplyr::select(gene, correlation.rnd, order2) %>% 
  dplyr::mutate(Status = "Random Correlation")

Spatcorrna2<- Spatcorrna2[order(Spatcorrna2$order2),] %>% 
  dplyr::rename(order1= order2) %>% 
  dplyr::rename(correlation = correlation.rnd)


Spatcorrnatot <- rbind(Spatcorrna1, Spatcorrna2)


ggplot(Spatcorrnatot, aes (x= order1, y= correlation, color = Status)) +
  geom_point()+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation')+
  xlab("Order")+
  ylab("Correlation")



##bulk RNA expression vs spatial RNA expression
bulkspatRNA <- readcount.vstDEGRADE %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate(gene_id = NULL) %>% 
  select(gene_name, everything()) %>% 
  inner_join(NanoString_RNA_data_col %>% dplyr::select(Gene), by = c("gene_name" = "Gene")) %>% 
  

bulkspatRNA <- bulkspatRNA[order(bulkspatRNA$gene_name), ]

spatRNAfilteredbulk <- bulkspatRNA %>% 
  dplyr::select(gene_name) %>% 
  left_join(NanoString_RNA_data_col, by = c("gene_name" = "Gene"))


bulkspatRNA2 <- bulkspatRNA %>% 
  as.data.frame(row.names = 1:nrow(.)) %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::mutate(averagebulk = rowMeans(.)) %>% 
  rownames_to_column("gene_name")

spatRNAfilteredbulk2 <- spatRNAfilteredbulk %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::mutate(averagespat = rowMeans(.)) %>% 
  rownames_to_column("gene_name")
  

exprbulkspat <- bulkspatRNA2 %>% 
  dplyr::select(gene_name, averagebulk) %>% 
  left_join(spatRNAfilteredbulk2 %>% dplyr::select(gene_name, averagespat))

rm(bulkspatRNA, spatRNAfilteredbulk, bulkspatRNA2, spatRNAfilteredbulk2)


ggplot(exprbulkspat, aes(averagebulk, log(averagespat)))+
  geom_point()+
  theme_bw()


cor.test(exprbulkspat$averagebulk, exprbulkspat$averagespat)

#DIFFERENTIAL EXPRESSION --> no use

NanoString_RNA_metadata <- NanoString_RNA_metadata %>% 
  dplyr::filter(Segment_tags == "Tumour_high")

DiffRNA_NanoString <- NanoString_RNA_data_col %>% 
  column_to_rownames("Gene") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id") %>% 
  dplyr::filter(patient_id %in% NanoString_RNA_metadata$pid_roi) %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame()


designRNANano <- model.matrix(~NanoString_RNA_metadata$Resection)
head(designRNANano)


fitRNANano <- lmFit(DiffRNA_NanoString, designRNANano)
fitRNANano <- eBayes(fitRNANano)
difRNANano <- topTable(fitRNANano,sort="none",n=Inf) %>% 
  rownames_to_column("Gene")




NanoString_protein_metadata <- NanoString_protein_metadata %>% 
  dplyr::filter(Segment_tags == "Tumour_high")

Diffprot_NanoString <- NanoString_protein_data_col %>% 
  column_to_rownames("Protein") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("patient_id") %>% 
  dplyr::filter(patient_id %in% NanoString_protein_metadata$pid_roi) %>% 
  column_to_rownames("patient_id") %>% 
  t() %>% 
  as.data.frame()




designprotNano <- model.matrix(~NanoString_protein_metadata$Resection)
head(designprotNano)


fitprotNano <- lmFit(Diffprot_NanoString, designprotNano)
fitprotNano <- eBayes(fitprotNano)
difprotNano <- topTable(fitprotNano,sort="none",n=Inf) %>% 
  rownames_to_column("Protein")



##OVERLAP DE BULK AND SPATIAL
spatbulkRNA<- deseq2ResultsGRADE %>% 
  inner_join(difRNANano, by = c("gene_name" = "Gene")) %>% 
  dplyr::mutate(group = "Not Significant") %>% 
  dplyr::mutate(group = ifelse(padj < 0.05, "Significant bulk RNA", group)) %>% 
  dplyr::mutate(group = ifelse(adj.P.Val < 0.05, "Significant spatial RNA", group)) %>% 
  dplyr::mutate(group = ifelse(padj < 0.05 & adj.P.Val < 0.05, "Significant both", group))








#1673 overlappen er 1494


ggplot(spatbulkRNA, aes(x = stat, y = t, color  = group, label = gene_name))+
  geom_point()+
  theme_bw()+
  xlab("bulk RNA")+
  ylab ("spatial RNA")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")

cor.test(spatbulkRNA$stat, spatbulkRNA$t)



#PROTEIN


spatbulkprot<- diffprot3 %>% 
  inner_join(difprotNano, by = c("Protein" = "Protein")) %>% 
  dplyr::mutate(group = "Not Significant") %>% 
  dplyr::mutate(group = ifelse(adj.P.Val.x < 0.05, "Significant bulk protein", group)) %>% 
  dplyr::mutate(group = ifelse(adj.P.Val.y < 0.05, "Significant spatial protein", group)) %>% 
  dplyr::mutate(group = ifelse(adj.P.Val.x < 0.05 & adj.P.Val.y < 0.05, "Significant protein", group))





ggplot(spatbulkprot, aes(x = t.x, y = t.y, color  = group, label = Protein))+
  geom_point()+
  theme_bw()+
  xlab("bulk protein")+
  ylab ("spatial protein")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")

cor.test(spatbulkprot$t.x, spatbulkprot$t.y)
