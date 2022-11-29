source("~/projects/glass/load_proteomics_data.R")
library(reshape2)
library(ggrepel)
library(fgsea)

dge.partially.paired.clusters <-readRDS("~/projects/glass/dge.partially.paired.clusters.Rds")
#highprotexpr2 has zero values, corproteinordered has imputed values
#totale correlatie RNA en protein----

##order RNA and protein data the same ----
sortpatientid <- sort(intersect(colnames(total_RNA), colnames(protein_counts_filtered)))

RNAordered <- total_RNA %>% 
  dplyr::select(c('gene_name', sortpatientid))

RNAordered <- RNAordered[order(RNAordered$gene_name),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL) 

rm(total_RNA)

proteinordered <- protein_counts_filtered %>% 
  dplyr::select(c('Protein', sortpatientid))

proteinordered <- proteinordered[order(proteinordered$Protein),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL) 

rm(protein_counts_filtered)

RNAcountsordered <- readcounts_filteredtotal %>% 
  dplyr::select(c('gene_name', sortpatientid))

RNAcountsordered <- RNAcountsordered[order(RNAcountsordered$gene_name),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL) 

rm(readcounts_filteredtotal)

##correlatieRNAProtein----
corproteinordered <- proteinordered %>% 
  column_to_rownames('Protein') %>% 
  t() %>% 
  as.data.frame()

rm(proteinordered)

corRNAordered <- RNAordered %>% 
  column_to_rownames('gene_name') %>% 
  t() %>% 
  as.data.frame() 

rm(RNAordered)

save(corRNAordered, file = "corRNAordered.Rdata")
save(RNAcountsordered, file = "RNAcountsordered.Rdata")
save(corproteinorderedNA, file = "corproteinorderedNA.Rdata")
save(corproteinorderedNATR, file = "corproteinorderedNATR.Rdata")

#high RNA expressed genes compared to high correlation

##correlatie RNA and protein compared to random -----
proteinorderedNA <- filteredNArawprotein %>% 
  rownames_to_column("Protein") %>% 
  dplyr::select(c('Protein', sortpatientid))

proteinorderedNA <- proteinorderedNA[order(proteinorderedNA$Protein),] %>% 
  tibble::rownames_to_column('RMMA') %>% 
  dplyr::mutate(RMMA=NULL) 

corproteinorderedNA <- proteinorderedNA %>% 
  column_to_rownames('Protein') %>% 
  t() %>% 
  as.data.frame() 

rm(proteinorderedNA, filteredNArawprotein, GLASS_filtered)

filteredentries_rawprotein <- corproteinorderedNA %>% 
  dplyr::mutate_all(function(x){return(ifelse(is.na(x),0,1))})


corproteinorderedNATR <- corproteinordered * filteredentries_rawprotein %>% 
  dplyr::mutate_all(function(x){return(ifelse(x == 0 ,NA_integer_,x))})

rm(filteredentries_rawprotein)

corRNAprotOLDSPEAR=data.frame()
for(i in colnames(corRNAordered)[1:length(corRNAordered)]) {
  print(i)
  a <- corRNAordered %>% dplyr::select(all_of(i))
  b <- corproteinordered %>% dplyr::select(all_of(i))
  
  
  print(paste0("Gene ",i,": ",cor(a,b)))
  
  corRNAprotOLDSPEAR <- rbind(corRNAprotOLDSPEAR, data.frame(gene=i, 
                                                             correlation = as.numeric(cor(a,b, method = "spearman"))
                                                            
  ))
}
rm(i,a,b)

corRNAprotOLDSPEAR <- corRNAprotOLDSPEAR[order(-corRNAprotOLDSPEAR$correlation),]


#spearman correlation transformed protein with NA values
corRNAprotNASPEAR=data.frame()
for(i in colnames(corRNAordered)[1:length(corRNAordered)]) {
  a <- corRNAordered %>% dplyr::select(all_of(i))
  b <- corproteinorderedNATR %>% dplyr::select(all_of(i))
  
  corRNAprotNASPEAR <- rbind(corRNAprotNASPEAR, data.frame(gene=i, 
                                                 correlation = as.numeric(cor(a,b, method = "spearman", use = "pairwise.complete.obs")),
                                                  correlation.rnd = as.numeric(cor(permute::shuffle(a),b, method = "spearman", use = "pairwise.complete.obs"))))
                                                 
                                              
}

corRNAprotNASPEAR <- corRNAprotNASPEAR[order(-corRNAprotNASPEAR$correlation),]

rm(corproteinorderedNA, i, a, b)

#mean and median of corRNAprotOLD and corRNAprotNA
#meancorRNAprotOLD <-mean(corRNAprotOLD$correlation)
#meancorRNAprotNA <- mean(corRNAprotNA$correlation)
#mediancorRNAprotOLD <- median(corRNAprotOLD$correlation)
#mediancorRNAprotNA <- median(corRNAprotNA$correlation)


corRNAprotNASPEAR %>% dplyr::filter(correlation > 0.7)
corRNAprotNASPEAR %>% dplyr::filter(correlation < 0.7 & correlation >0.5)
corRNAprotNASPEAR %>% dplyr::filter(correlation < 0.5 & correlation > 0.3)
corRNAprotNASPEAR %>% dplyr::filter(correlation <0.3 & correlation >= 0 )

corRNAprotNASPEAR %>% dplyr::filter(correlation < -0.5)
corRNAprotNASPEAR %>% dplyr::filter(correlation >-0.5 & correlation < -0.3)
corRNAprotNASPEAR %>% dplyr::filter(correlation >-0.3 & correlation <0)

#genen met de hoogste correlatie corplotten
corRNAprothigh <- corRNAprotNASPEAR %>% 
  dplyr::mutate(correlation.rnd=NULL) %>% 
  dplyr::select(gene)

rm(a,b,i, removed, a2, b2)

corrndcompare <- corRNAprotNASPEAR %>%
  dplyr::mutate(order1 = rank(-correlation, ties.method = "first")) %>% 
  dplyr::mutate(order2 = rank(-correlation.rnd, ties.method = "first")) 

corrna1 <- corrndcompare %>% 
  dplyr::select(gene, correlation,order1) %>% 
  dplyr::mutate(Status = "Correlation")

corrna2 <- corrndcompare %>% 
  dplyr::select(gene, correlation.rnd, order2) %>% 
  dplyr::mutate(Status = "Random Correlation")

corrna2<- corrna2[order(corrna2$order2),] %>% 
  dplyr::rename(order1= order2) %>% 
  dplyr::rename(correlation = correlation.rnd)

corrna3 <- corRNAprotOLDSPEAR %>% 
  dplyr::mutate(order1 = rank(-correlation, ties.method = "first")) %>% 
  dplyr::select(gene, correlation,order1) %>% 
  dplyr::mutate(Status = "Imputed Correlation")

corrnatot <- rbind(corrna1, corrna2, corrna3)


save (corrnatot, file = "corrnatot.Rdata")
# corrnatot2 <- tibble(corrnatot) %>%  dplyr::left_join(tibble(dge.partially.paired.clusters), by=c('gene'='gene_name')) %>% 
#   dplyr::filter(!is.na(gene_uid)) %>% 
#   dplyr::filter(Status == "Correlation")


# 
# #wilcoxon signed rank test
# Realcor <- corrna3$correlation
# Randomcor <- corrna2$correlation
# realvsrandom <- data.frame(
#   group = rep(c("Real Correlation", "Random Correlation"), each = 3135),
#   correlation = c(Realcor, Randomcor))
# 
# dplyr::group_by(realvsrandom, group) %>% 
#   summarise(
#     count = n(),
#     median = median(correlation),
#     IQR = IQR(correlation)
#   )
# 
# library("ggpubr")
# ggboxplot(realvsrandom, x = "group", y = "correlation", 
#           color = "group", palette = c("#00AFBB", "#E7B800"),
#           ylab = "Correlation", xlab = "Groups")
# 
# wilcoxrealrandom <- wilcox.test(Realcor, Randomcor)
# 
# rm(wilcoxrealrandom, realvsrandom)
# 

# #GROUPING
# Yhighrnaexprclust <- as.data.frame(intersect(names(corproteinordered) , dge.partially.paired.clusters$gene_name)) %>% 
#   dplyr::rename(Overlap = "intersect(names(corproteinordered), dge.partially.paired.clusters$gene_name)")
# 
# corrnatot$group <- "0"
# corrnatot$group[corrnatot$gene %in% c("LMNB1", "ACTC1", "CDK1", "H4C1", "H2AC4", "H2AC11", "H2BC12", "H1-5", "H2BC17")] <- "up.1"
# corrnatot$group[corrnatot$gene %in% c("COL6A3", "COL6A2", "TGFBI","HSPG2","LUM","FCGBP","CHI3L1","COL1A2","ACAN","VGF","COL4A2","COL1A1")] <- "up.2"
# #corrnatot$group[corrnatot$gene %in% c("AQP1", "FABP5")] <- "up.3"
# corrnatot$group[corrnatot$gene %in% c("MPP6", "OGFRL1", "GJA1", "RGS6")] <- "down"
# 
# corrnatot <- corrnatot %>%  dplyr::mutate(colors = ifelse(group == "0", Status, group)) %>% 
#   dplyr::mutate(colors = ifelse(Status == "Random Correlation" & colors != "Random Correlation", Status, colors)) %>% 
#   dplyr::mutate(colors = ifelse(Status == "Correlation" & colors != "Correlation", Status, colors))
# 
# ##PAPER
# 
# corrnatotPAPER <- rbind(corrna1, corrna2)
# corrnatotPAPER$group <- "0"
# corrnatotPAPER$group[corrnatotPAPER$gene %in% c("LMNB1", "ACTC1", "CDK1", "H4C1", "H2AC4", "H2AC11", "H2BC12", "H1-5", "H2BC17")] <- "up.1"
# corrnatotPAPER$group[corrnatotPAPER$gene %in% c("COL6A3", "COL6A2", "TGFBI","HSPG2","LUM","FCGBP","CHI3L1","COL1A2","ACAN","VGF","COL4A2","COL1A1")] <- "up.2"
# #corrnatot$group[corrnatot$gene %in% c("AQP1", "FABP5")] <- "up.3"
# corrnatotPAPER$group[corrnatotPAPER$gene %in% c("MPP6", "OGFRL1", "GJA1", "RGS6")] <- "down"
# 
# corrnatotPAPER <- corrnatotPAPER %>%  dplyr::mutate(colors = ifelse(group == "0", Status, group)) %>% 
#   dplyr::mutate(colors = ifelse(Status == "Random Correlation" & colors != "Random Correlation", Status, colors))
# 
# save(corrnatotPAPER, file= "corrnatotPAPER.Rdata")
