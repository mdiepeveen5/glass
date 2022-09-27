#needs protupprot --> protein significant omhoog en RNA gematcht, maar RNA hoeft niet significant te zijn.


#PROTEIN significant omhoog, hangen aan RNA expressie en kijken of je ook cell cycle en periciet clusters vindt

load(file = "readcount.vstDE.Rdata")

RNAexprprotsignup <-  readcount.vstDE %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate(gene_id=NULL) %>% 
  select(gene_name, everything()) %>% 
  right_join(diffprot2signup %>% dplyr::select(Protein), by = c("gene_name" = "Protein")) %>% 
  column_to_rownames("gene_name")
  
hclustprot <- hclust(as.dist(1-cor(t(RNAexprprotsignup2))))
o <- hclustprot$labels

corrplot(cor(as.data.frame(t(RNAexprprotsignup2))), order = "hclust", tl.cex = 0.4)
plot(hclustprot)

clustersprot1 <-cutree(hclustprot, k = 12) %>% 
  purrr::keep(function(x) x== "1") %>% names

print(as.data.frame(clustersprot1), row.names = FALSE)

clustersprot2 <-cutree(hclustprot, k = 12) %>% 
  purrr::keep(function(x) x== "5") %>% names
print(as.data.frame(clustersprot2), row.names = FALSE)



RNAexprprotsignup2 <-RNAexprprotsignup %>% rownames_to_column("gene_id") %>% 
  rbind(readcount.vstDE %>% rownames_to_column("gene_id") %>% dplyr::filter(gene_id == "ENSG00000174807.4")) %>% 
  dplyr::mutate(gene_id = ifelse (gene_id == "ENSG00000174807.4", "CD248", gene_id)) %>% 
  column_to_rownames("gene_id")


cellcyclecollagen <- c("COL6A3", "TGFBI", "TXNDC5", "SERPINH1",
                       "BAX","COL18A1","COL6A2","BGN","RCC2","RCC1","H1-5","MCM3","MCM7","PRPF4","SSRP1","JPT1", "CD248") %>% 
  as.data.frame() %>% 
  rename("." = "gene_name")




collagenandcellcyle <- RNAexprprotsignup2 %>%
  rownames_to_column("gene_name") %>% 
  dplyr::right_join(cellcyclecollagen) %>% 
  distinct() %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("gene_name")

corrplot(cor(as.data.frame(t(collagenandcellcyle))), order = "hclust", tl.cex = 0.4)


RNAexprprotsignup5<- RNAexprprotsignup%>% 
  rownames_to_column("gene_name")

print(as.data.frame(RNAexprprotsignup5$gene_name), row.names = F)

##COLLAGEN : https://biit.cs.ut.ee/gplink/l/WfPQK4q1QG


load(file = "Protallpatients.Rdata")

##CELL CYCLE : 

corplotdiffprotsignup <- Protallpatients %>% 
  right_join(diffprot2signup %>% dplyr::select(Protein)) %>% 
  column_to_rownames("Protein")
  
hclustprot <- hclust(as.dist(1-cor(t(corplotdiffprotsignup))))
corrplot(cor(as.data.frame(t(corplotdiffprotsignup))), order = "hclust", tl.cex = 0.4, type = "lower")
plot(hclustprot)


cyclecolup1 <- c("COL6A3", "TGFBI", "TXNDC5", "SERPINH1",
  "BAX","COL18A1","COL6A2","BGN","RCC2","RCC1","H1-5","MCM3",
  "MCM7","PRPF4","SSRP1","JPT1", "CD248", "LMNB1", "H4C1","H2AC4", "H2AC11",
  "H2BC12", "H2BC17", "CDK1", "ACTC1", "HSPG2", "CHI3L1", "COL1A2", 
  "VGF","LUM", "COL4A2", "ACAN", "COL1A1", "FCGBP") %>% 
  as.data.frame() %>% 
  rename("." = "gene_name")

cellcycleprotup1 <- c("PCNA", "MCM7", "MDC1", "RPA3", "ASNS", "RPRD1B","CDC73", "SOX2", "MYBBP1A",
                      "SMARCA5", "GOLGA2", "RAD21", "MCM3", "ZNF207", "PDS5A", "PSME3", "ACTL6A",
                      "MAPK14", "RCC2","LIG3","MRE11", "RBBP4", "RCC1","NUP43","PRPF40A","HCFC1",
                      "BAX", "EEF1AKNMT", "LMNB1", "H4C1","H2AC4", "H2AC11",
                      "H2BC12", "H2BC17", "CDK1", "ACTC1", "H1-5")%>% 
  as.data.frame() %>% 
  rename("." = "gene_name")


#####################################################

DEcorplRNA <-  readcount.vstDE %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gene.annot2) %>% 
  dplyr::mutate(gene_id=NULL) %>% 
  select(gene_name, everything()) %>% 
  right_join(cellcycleprotup1) %>% 
  column_to_rownames("gene_name")


load(file = "Protallpatients.Rdata")

DEcorplprot <- Protallpatients %>% 
  right_join(cellcycleprotup1, by = c("Protein" = "gene_name")) %>% 
  column_to_rownames("Protein")


hclustcorpl <- hclust(as.dist(1-cor(t(DEcorplprot))))
ordercorpl <- hclustcorpl$labels[hclustcorpl$order]


corprot <- cor(as.data.frame(t(DEcorplprot)) %>%  dplyr::select(ordercorpl), use = "pairwise.complete.obs")

corrplot(corprot)

corRNA <- cor(as.data.frame(t(DEcorplRNA)) %>%  dplyr::select(ordercorpl), use = "pairwise.complete.obs")

corprot[upper.tri(corprot)] <- corRNA[upper.tri(corprot)]

labelsprotclust = c("red","black","black","red","red","black","black","black","red",  
                    "black","black","black","black","black","black","black","black","black",  
                    "black","black","black","black","black","black","black","black","black",   
                    "black","red","red","red","red","red","black","black","black",  
                    "black") 


corrplot(corprot, order= 'original',tl.cex = 0.4, diag = F, tl.col = labelsprotclust)



hclustcorpl2 <- hclust(as.dist(1-cor(t(DEcorplRNA))))
ordercorpl2 <- hclustcorpl2$labels[hclustcorpl2$order]


corprot2 <- cor(as.data.frame(t(DEcorplprot)) %>%  dplyr::select(ordercorpl2), use = "pairwise.complete.obs")

corrplot(corprot2)

corRNA2 <- cor(as.data.frame(t(DEcorplRNA)) %>%  dplyr::select(ordercorpl2), use = "pairwise.complete.obs")

corprot2[upper.tri(corprot2)] <- corRNA2[upper.tri(corprot2)]

labelsRNAclust = c("black","black","black","black","black","black","black","black","black",  
                    "black","black","black","black","black","black","black","black","black",  
                    "red","black","red","red","red","red","red","red","red",   
                    "red","black","black","black","black","black","black","black","black",  
                    "black") 


corrplot(corprot2, order= 'original',tl.cex = 0.4, diag = F, tl.col = labelsRNAclust)



rm(hclustcorpl, hclustcorpl2, labelsRNAclust, labelsprotclust, ordercorpl, ordercorpl2, DEcorplprot, DEcorplRNA, corRNA,
   corRNA2, corprot, corprot2)


