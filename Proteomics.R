load(file = "corrnatot.Rdata")


library("ggpmisc")
Yhighrnaexprclust <- as.data.frame(intersect(names(corproteinordered) , dge.partially.paired.clusters$gene_name)) %>% 
  dplyr::rename(Overlap = "intersect(names(corproteinordered), dge.partially.paired.clusters$gene_name)")

corrnatot$group <- "0"
#corrnatot$group[corrnatot$gene %in% c("LMNB1", "ACTC1", "CDK1", "H4C1", "H2AC4", "H2AC11", "H2BC12", "H1-5", "H2BC17")] <- "up.1"
#corrnatot$group[corrnatot$gene %in% c("COL6A3", "COL6A2", "TGFBI","HSPG2","LUM","FCGBP","CHI3L1","COL1A2","ACAN","VGF","COL4A2","COL1A1")] <- "up.2"
#corrnatot$group[corrnatot$gene %in% c("AQP1", "FABP5")] <- "up.3"
#corrnatot$group[corrnatot$gene %in% c("MPP6", "OGFRL1", "GJA1", "RGS6")] <- "down"

corrnatot <- corrnatot %>%  dplyr::mutate(colors = ifelse(group == "0", Status, group)) %>% 
  dplyr::mutate(colors = ifelse(Status == "Random Correlation" & colors != "Random Correlation", Status, colors)) %>% 
  dplyr::mutate(colors = ifelse(Status == "Correlation" & colors != "Correlation", Status, colors))

#plot correlation, random correlation and imputed correlation
ggplot(corrnatot, aes (x= order1, y= correlation, color = colors, label = gene)) +
  scale_color_manual(values=c("chartreuse","blue4", "yellow")) + 
  geom_point(data = subset(corrnatot, group == '0'), size=2) +
  geom_point(data = subset(corrnatot, group != '0'), size=2) +
  geom_text_repel(data = subset(corrnatot, group != '0' & Status != 'Correlation' & Status != "Random Correlation" )) +
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")

rm(corrna1, corrna2, corrna3, corrndcompare)


ggplot(corrnatot, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(corrnatot, group == '0' & Status != "Imputed Correlation"), size=2) +
  scale_color_manual(values=c("chartreuse","blue4", "yellow"))+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation')+
  xlab("Order")+
  ylab("Correlation")


#Example high, low and anticorrelation
#High correlation
highcorTNR = cbind(as.data.frame(corRNAordered$TNR), as.data.frame(corproteinorderedNATR$TNR))
ggplot(highcorGSTM3, aes(x = corRNAordered$TNR, y = corproteinorderedNATR$TNR) )+
  geom_point() +
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title="Correlation of TNR",
        x ="RNA expression", y = "Protein expression")

#No correlation
nocorFSD1 = cbind(as.data.frame(corRNAordered$FSD1), as.data.frame(corproteinorderedNATR$FSD1))
ggplot(nocorFSD1, aes(x = corRNAordered$FSD1, y = corproteinorderedNATR$FSD1) )+
  geom_point() +
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title="Correlation of FSD1",
        x ="RNA expression", y = "Protein expression")

anticorNDUFA13 = cbind(as.data.frame(corRNAordered$NDUFA13), as.data.frame(corproteinorderedNATR$NDUFA13))
ggplot(anticorNDUFA13, aes(x = corRNAordered$NDUFA13, y = corproteinorderedNATR$NDUFA13) )+
  geom_point() +
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title="Correlation of NDUFA13",
        x ="RNA expression", y = "Protein expression")


#RNA correlation with protein and protein with that RNA
corprotCTNND2RNAHLAA <- cbind(as.data.frame(corRNAordered$`HLA-A`), as.data.frame(corproteinorderedNATR$CTNND2))
ggplot(corprotCTNND2RNAHLAA, aes(x = corRNAordered$`HLA-A`, y = corproteinorderedNATR$CTNND2) )+
  geom_point() +
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title="Correlation of RNA HLA-A and protein CTNND2",
        x ="RNA expression", y = "Protein expression")

corprotHLAARNACTNND2 <- cbind(as.data.frame(corRNAordered$CTNND2), as.data.frame(corproteinorderedNATR$`HLA-A`))
ggplot(corprotCTNND2RNAHLAA, aes(x = corRNAordered$CTNND2, y = corproteinorderedNATR$`HLA-A`) )+
  geom_point() +
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title="Correlation of RNA CTNND2 and protein HLA-A",
        x ="RNA expression", y = "Protein expression")

#you can load the datasets from savings
##cor RNA expressie----
highRNAexpr <- t(corRNAordered) %>% 
  as.data.frame() %>% 
  dplyr::mutate(sumrow = rowSums(.))
highRNAexpr <- highRNAexpr[order(-highRNAexpr$sumrow),] 

highRNAexprplot <- highRNAexpr %>% 
  dplyr::select(sumrow) %>% 
  tibble::rownames_to_column("gene_name")

highRNAcor <- left_join(highRNAexprplot, corRNAprotNASPEAR, by = c("gene_name" = "gene")) %>% 
  dplyr::mutate(correlation.rnd = NULL)

ggplot(data = highRNAcor, aes (x = sumrow, y = correlation))+
  geom_point()+
  labs (title = "Correlation increase with RNA expression", x= "Sum of RNA expression per gene", y = "Correlation")+
  stat_poly_eq()+
  stat_poly_line()+
  theme_bw()

cor.test(log( 1+highRNAcor$sumrow), highRNAcor$correlation)

rm(highRNAexprplot, highRNAcor, highRNAexpr)



##sumrawcountplots----
RNAcountsordered <- RNAcountsordered %>% 
  column_to_rownames('gene_name') %>% 
  t() %>% 
  as.data.frame() 

highrawRNAcounts <- t(RNAcountsordered) %>% 
  as.data.frame() %>% 
  dplyr::mutate(meanrow = rowMeans(.))

highrawRNAcounts <- highrawRNAcounts[order(-highrawRNAcounts$meanrow),] 

highRNAcountsplot <- highrawRNAcounts %>% 
  dplyr::select(meanrow) %>% 
  tibble::rownames_to_column("gene_name")

highRNAcountcor <- left_join(highRNAcountsplot, corRNAprotNASPEAR, by = c("gene_name" = "gene")) %>% 
  dplyr::mutate(correlation.rnd = NULL)

ggplot(data = highRNAcountcor, aes (x = log(1+meanrow), y = correlation))+
  geom_point()+
  geom_smooth(method = lm) +
  labs ()+
  xlab ("Meanrow of log (1+ raw RNA count)) expression per gene")+
  ylab ("Correlation Protein x RNA")

cor.test(log( 1+highRNAcountcor$meanrow), highRNAcountcor$correlation)

rm(highRNAcountsplot, highRNAcountcor, highrawRNAcounts)


##sum zero count plots ----

highrawRNAzerocounts <- t(RNAcountsordered) %>% 
  as.data.frame()

highrawRNAzerocounts$zeros <- rowSums(highrawRNAzerocounts == 0)

highrawRNAzerocounts <- highrawRNAzerocounts[order(-highrawRNAzerocounts$zeros),] 

highRNAzerocountsplot <- highrawRNAzerocounts %>% 
  dplyr::select(zeros) %>% 
  tibble::rownames_to_column("gene_name")

highRNAzerocountcor <- left_join(highRNAzerocountsplot, corRNAprotNASPEAR, by = c("gene_name" = "gene")) %>% 
  dplyr::mutate(correlation.rnd = NULL)

ggplot(data = highRNAzerocountcor, aes (x = log(1+zeros), y = correlation))+
  geom_point()+
  labs (title = "Correlation increase with RNA zero counts", x= "Number of zero counts per gene (1+log(counts))", y = "Correlation")+
  stat_poly_eq()+
  stat_poly_line()+
  theme_bw()

cor.test(log( 1+highRNAzerocountcor$zeros), highRNAzerocountcor$correlation)


rm(highRNAzerocountcor, highRNAzerocountsplot, highrawRNAzerocounts)

##cor protein expressie----

filteredentries_rawprotein <- GLASS_filtered %>% 
  dplyr::mutate_all(function(x){return(ifelse(x == "Filtered",NA,1))})
#highprotexpr with NA = 0
highprotexpr <- ((t(corproteinorderedNATR))) %>% 
  as.data.frame() %>% 
  dplyr::mutate(meanrow = rowMeans(.))
highprotexpr  <-highprotexpr[order(-highprotexpr$meanrow),] 

highprotexprplot <- highprotexpr %>% 
  dplyr::select(meanrow) %>% 
  tibble::rownames_to_column("protein_name")

highprotcor <- left_join(highprotexprplot, corRNAprotNASPEAR, by = c("protein_name" = "gene")) %>% 
  dplyr::mutate(correlation.rnd = NULL)

ggplot(data = highprotcor, aes (x = meanrow, y = correlation))+
  geom_point()+
  labs (title = "Correlation increase with protein expression", x= "Mean protein expression per gene", y = "Correlation")+
  stat_poly_eq()+
  stat_poly_line()+
  theme_bw()


cor.test(highprotcor$meanrow, highprotcor$correlation)
rm(highprotexprplot, highprotcor, GLASS_filtered)

#zerocounts summation protein
protzerocounts <- ((t(corproteinorderedNATR))) %>% 
  as.data.frame() %>% 
  dplyr::mutate_all(function(x){return (ifelse(is.na(x), 0,x))}) %>% 
  dplyr::mutate_all(function(x) as.numeric(as.character(x)))

protzerocounts $zeros <- rowSums(protzerocounts == 0)
protzerocounts  <-protzerocounts[order(-protzerocounts$zeros),] 

protzerocountsplot <- protzerocounts %>% 
  dplyr::select(zeros) %>% 
  tibble::rownames_to_column("gene_name")

protzerocountscor <- left_join(protzerocountsplot, corRNAprotNASPEAR, by = c("gene_name" = "gene")) %>% 
  dplyr::mutate(correlation.rnd = NULL)

ggplot(data = protzerocountscor, aes (x = log(1+zeros), y = correlation))+
  geom_point()+
  labs (title = "Correlation increase with protein zero counts", x= "Number of zero counts per gene (1+log(counts))", y = "Correlation")+
  stat_poly_eq()+
  stat_poly_line()+
  theme_bw()

cor.test(protzerocountscor$zeros, protzerocountscor$correlation)
rm(protzerocounts, protzerocountsplot, filteredentries_rawprotein, protzerocountscor)

#correlatieplot op genen met hoogste expressie
corplotRNAhigh <- t(highRNAexpr) %>% 
  as.data.frame() %>% 
  dplyr::mutate(sumrow=NULL) %>% 
  slice_head(n=20)

corrplot::corrplot(cor(corplotRNAhigh[1:20,1:20]))

corplotprothigh <- t(highprotexpr) %>% 
  as.data.frame() %>% 
  dplyr::mutate(meanrow=NULL) %>% 
  slice_head(n=20)

corrplot::corrplot(cor(corplotprothigh[1:20,1:20]))

rm(corplotprothigh, corplotRNAhigh)

#make meanrow and sumrow NULL
highRNAexpr2 <- highRNAexpr %>% 
  dplyr::mutate(sumrow=NULL)
highprotexpr2 <- highprotexpr %>% 
  dplyr::mutate(meanrow=NULL)

#h <- hclust(as.dist(1-cor(t(highRNAexpr2[1:100,]))))
#o <- h$labels
corrplot(cor(as.data.frame(t(highRNAexpr2[1:100,])) %>%  dplyr::select(o)) , order="hclust")

corrplot(cor(as.data.frame(t(highRNAexpr2[1:100,]))), order = "hclust")

corrplot(cor(as.data.frame(t(highprotexpr2[1:100,]))), type = "upper")

#protein highest cor
protRNAprothigh <- t(corproteinordered) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("protein")

protRNAprothigh <- dplyr::inner_join(corRNAprothigh,protRNAprothigh, by = c("gene" = "protein")) 

rownames(protRNAprothigh) <- NULL

protRNAprothigh<- protRNAprothigh %>% 
  tibble::column_to_rownames("gene") 

#NON IMPUTED VALUES
protRNAprothighNA <- t(corproteinorderedNATR) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("protein")

protRNAprothighNA <- dplyr::inner_join(corRNAprothigh,protRNAprothighNA, by = c("gene" = "protein"))  

rownames(protRNAprothighNA) <- NULL

protRNAprothighNA<- protRNAprothighNA %>% 
  tibble::column_to_rownames("gene") 

#RNA highest cor
RNARNAprothigh <- t(corRNAordered)%>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene")

 
RNARNAprothigh <- dplyr::inner_join(corRNAprothigh, RNARNAprothigh, by = c("gene" = "gene"))

rownames(RNARNAprothigh) <- NULL

RNARNAprothigh <- RNARNAprothigh %>% 
  tibble::column_to_rownames("gene")


#corrplot RNA protein
corrplot(cor(as.data.frame(t(protRNAprothighNA[1:100,])), as.data.frame(t(RNARNAprothigh[1:100,])), 
             use = "pairwise.complete.obs", method = "spearman"), order = "hclust", tl.cex = 0.4)
mtext(text = "Proteome", side = 2, line = -19, cex = 3.6)
mtext(text = "Transcriptome", side = 1, line =-47.5, cex = 3.6)


#corrplot RNARNA
corrplot(cor(as.data.frame(t(RNARNAprothigh[1:100,])), use = "pairwise.complete.obs", method = "spearman"), order = "hclust", tl.cex = 0.4)
RNAcor <-(cor(as.data.frame(t(RNARNAprothigh)), use = "pairwise.complete.obs", method = "spearman"))
mean(cor(as.data.frame(t(RNARNAprothigh)), use = "pairwise.complete.obs", method = "spearman"))
median(cor(as.data.frame(t(RNARNAprothigh)), use = "pairwise.complete.obs", method = "spearman"))

#corrplot proteinprotein
corrplot(cor(as.data.frame(t(protRNAprothighNA[1:100,])), use = "pairwise.complete.obs", method = "spearman"), order = "hclust", tl.cex = 0.4)
mean(cor(as.data.frame(t(protRNAprothighNA)), use = "pairwise.complete.obs", method = "spearman"))
median(cor(as.data.frame(t(protRNAprothighNA)), use = "pairwise.complete.obs", method = "spearman"))

rm(RNAcor)


#make triangular corrplot PROTEIN RNA
#hclustering 
#protRNA
hclustprotRNA <- hclust(as.dist(1-cor(as.data.frame(t(protRNAprothighNA[1:100,])), as.data.frame(t(RNARNAprothigh[1:100,])), 
                                      use = "pairwise.complete.obs", method = "spearman")))
orderprotRNA <- hclustprotRNA$labels[hclustprotRNA$order]

plot(hclustprotRNA)

#see if clusters are enriched in certain cell types --> answer is NO
clustersprotRNA1 <-cutree(hclustprotRNA, k = 3) %>% 
  purrr::keep(function(x) x== "1") %>% names

print(as.data.frame(clustersprotRNA1), row.names = FALSE)


clustersprotRNA2 <- cutree(hclustprotRNA, k=3) %>% 
  purrr::keep(function(x) x== "2") %>% names

print(as.data.frame(clustersprotRNA2), row.names = FALSE)

clustersprotRNA3 <- cutree(hclustprotRNA, k=3) %>% 
  purrr::keep(function(x) x== "3") %>% names

print(as.data.frame(clustersprotRNA3), row.names = FALSE)

#https://biit.cs.ut.ee/gplink/l/9PHS2hPtQX

print(as.data.frame(corRNAprotNASPEAR$gene[1:100]), row.names = FALSE)



clustersprotRNA3 <- dendextend::prune(hclustprotRNA,cutree(hclustprotRNA,3) %>% 
                                        purrr::keep(function(x) x== "3") %>% names)

#RNA
hclustRNA <- hclust(as.dist(1-cor(t(RNARNAprothigh[1:100,]))))
orderRNA <- hclustRNA$labels[hclustRNA$order]

#protein
hclustprot <- hclust(as.dist(1-cor(t(protRNAprothighNA[1:100,]), use = "pairwise.complete.obs")))
orderprot <- hclustprot$labels[hclustprot$order]



corrplot(cor(as.data.frame(t(protRNAprothighNA[1:100,])) %>%  dplyr::select(orderprotRNA), use = "pairwise.complete.obs") , order= 'original',  type = "lower",tl.cex = 0.4)
corrplot(cor(as.data.frame(t(RNARNAprothigh[1:100,])) %>%  dplyr::select(orderprotRNA), use = "pairwise.complete.obs") , order= 'original',  type = "upper",tl.cex = 0.4, add =T)

corprot <- cor(as.data.frame(t(protRNAprothighNA[1:100,])) %>%  dplyr::select(orderprotRNA), use = "pairwise.complete.obs")

corRNA <- cor(as.data.frame(t(RNARNAprothigh[1:100,])) %>%  dplyr::select(orderprotRNA), use = "pairwise.complete.obs")

corprot[upper.tri(corprot)] <- corRNA[upper.tri(corprot)]

corrplot(corprot, order= 'original',tl.cex = 0.4, diag = F)


corrplot(cor(as.data.frame(t(protRNAprothighNA[1:100,])), use = "pairwise.complete.obs"), type = "upper", diag = F, tl.cex = 0.4)

corrplot(cor(as.data.frame(t(RNARNAprothigh[1:100,])), use = "pairwise.complete.obs"), type = "lower", diag = F, tl.cex = 0.4, add = T)


rm(protRNAprothigh, RNARNAprothigh, highprotexpr, highRNAexpr, corRNAprotOLD)


rm(orderRNA, orderprot, orderprotRNA, corRNA, corprot)

rm(hclustprot, hclustprotRNA, hclustRNA)



#dplyr::inner_join for combined RNA and protein data
#most abundant proteins selected on MAD
mostabundant_prot<- tmpprot %>% 
  tibble::column_to_rownames("Protein") %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=20) %>% 
  dplyr::mutate(mad=NULL)

