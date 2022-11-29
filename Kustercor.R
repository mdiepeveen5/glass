library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

data_koster  <- read_excel_allsheets("Table_EV1.xlsx")
data_koster2 <- read_excel_allsheets("Table_EV2.xlsx")


genes_koster <- as.data.frame(data_koster2$`B. Genes`) %>% 
  dplyr::select("Gene name", "Brain") %>% 
  inner_join(interprotRNA, by = c("Gene name" = "gene_name") )

proteins_koster <- as.data.frame(data_koster$`C. Genes`)%>% 
  dplyr::select("Gene name", "Brain") %>% 
  inner_join(interprotRNA, by = c("Gene name" = "gene_name") )

prot_gene_koster <- inner_join(genes_koster, proteins_koster, by = c("Gene name")) %>% 
  dplyr::filter(Brain.x != 0) %>% 
  dplyr::filter(Brain.y != 0) 

duplicated <- prot_gene_koster$`Gene name`[duplicated(prot_gene_koster$`Gene name`)]
  
prot_gene_koster <- prot_gene_koster %>% 
  dplyr::filter(`Gene name` %in% duplicated == F) %>% 
  dplyr::rename("Gene" = `Gene name`) %>% 
  inner_join(corRNAprotNASPEAR %>% dplyr::select(gene, correlation), by = c("Gene" = "gene")) %>% 
  dplyr::mutate(colours =ifelse(correlation > 0.5,  "High Correlation", "Low Correlation")) %>% 
  dplyr::mutate(colours = ifelse(correlation < -0.5, "Anti Correlation", colours))


ggplot(prot_gene_koster, aes(x = log(Brain.x), y = log(Brain.y), color = colours,  label = Gene))+
  geom_point(data = subset(prot_gene_koster, colours == "Low Correlation"))+
  geom_point(data = subset(prot_gene_koster, colours == "High Correlation"))+
  geom_point(data = subset(prot_gene_koster, colours == "Anti Correlation"))+
  theme_bw()+
  labs (title = "RNA x Protein", x= "RNA" , y = "Protein")+
  geom_text_repel(data = subset(prot_gene_koster, colours == "High Correlation"), color = "black")+
  scale_color_manual(values = c("grey", "red", "blue"),
                     name = "Correlation",
                     breaks=c("Low Correlation", "High Correlation", "Anti Correlation"),
                     labels = c("Low Correlation", "High Correlation", "Anti Correlation"))


sumRNAexpr <- corRNAordered %>% 
  t() %>% 
  as.data.frame %>% 
  dplyr::mutate(sumRNArow = rowSums(.)) %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sumRNArow)
  
  

sumprotexpr <- corproteinordered%>% 
  t() %>% 
  as.data.frame %>% 
  dplyr::mutate(sumprotrow = rowSums(.)) %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, sumprotrow) %>% 
  inner_join(sumRNAexpr)

ggplot(sumprotexpr, aes(x = sumRNArow, y = sumprotrow, label = gene))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title = "RNA x Protein", x= "RNA" , y = "Protein")+
  geom_text_repel()


##

genes_koster_total <- as.data.frame(data_koster2$`B. Genes`) %>% 
  dplyr::mutate(`Tissue enriched` = NULL) %>% 
  dplyr::mutate(`Group enriched` = NULL) %>% 
  dplyr::mutate(`Classification` = NULL) %>% 
  dplyr::mutate(`Gene ID` = NULL) %>% 
  dplyr::mutate(`Tissue enhance` = NULL) %>% 
  inner_join(interprotRNA, by = c("Gene name" = "gene_name")) 
  

proteins_koster_total <- as.data.frame(data_koster$`C. Genes`)%>% 
  dplyr::mutate(`Tissue enriched` = NULL) %>% 
  dplyr::mutate(`Group enriched` = NULL) %>% 
  dplyr::mutate(`Classification` = NULL) %>% 
  dplyr::mutate(`Gene ID` = NULL) %>% 
  dplyr::mutate(`Tissue enhanced` = NULL) %>% 
  inner_join(interprotRNA, by = c("Gene name" = "gene_name")) 

duplicated2 <- genes_koster_total$`Gene name`[duplicated(genes_koster_total$`Gene name`)]
duplicated3 <- proteins_koster_total$`Gene name`[duplicated(proteins_koster_total$`Gene name`)]

genes_koster_total <- genes_koster_total%>% 
  dplyr::filter(`Gene name` %in% duplicated2 == F) 


proteins_koster_total <- proteins_koster_total%>% 
  dplyr::filter(`Gene name` %in% duplicated3 == F)

#filter the datasets
genes_koster_total <- genes_koster_total %>% 
  dplyr::filter(`Gene name` %in% proteins_koster_total$`Gene name`)

proteins_koster_total <- proteins_koster_total%>% 
  dplyr::filter(`Gene name` %in% genes_koster_total$`Gene name`) 

#order datasets
genes_koster_total <- genes_koster_total[order(genes_koster_total$`Gene name`),] 

proteins_koster_total <- proteins_koster_total[order(proteins_koster_total$`Gene name`),] 
  
#make rownames zero so you can tibble the data 
row.names(genes_koster_total) <- NULL
row.names(proteins_koster_total) <- NULL

genes_koster_total <- genes_koster_total %>% 
  column_to_rownames("Gene name") %>% 
  t() %>% 
  as.data.frame()

proteins_koster_total <- proteins_koster_total %>% 
  column_to_rownames("Gene name") %>% 
  t() %>% 
  as.data.frame()


corKUSTER=data.frame()
for(i in colnames(genes_koster_total)[1:length(genes_koster_total)]) {
  a <- genes_koster_total %>% dplyr::select(all_of(i))
  b <- proteins_koster_total %>% dplyr::select(all_of(i))
  
  corKUSTER <- rbind(corKUSTER, data.frame(gene=i, 
                                                    correlation = as.numeric(cor(a,b, method = "spearman", use = "pairwise.complete.obs"))
  ))

  
}

corKUSTER <- corKUSTER[order(-corKUSTER$correlation),]


KUSTERcompare <- corKUSTER %>%
  dplyr::mutate(order1 = rank(-correlation, ties.method = "first")) %>% 
  
  dplyr::mutate(colours = ifelse(gene %in% corRNAprothigh$gene[1:100], "High", "Low"))



ggplot(KUSTERcompare, aes (x= order1, y= correlation, color = colours, label = gene)) +
  geom_point(data = subset(KUSTERcompare, colours == "Low"), size = 1)+
  geom_point(data = subset(KUSTERcompare, colours == "High"), size = 2)+
  geom_text_repel(data = subset(KUSTERcompare, colours == "High", max.overlaps = 10))+
  theme_bw()+
  scale_color_manual(values=c("red3","blue4")) +
  labs(title = "Correlation between RNA and protein among 27 measured tissues", color='Correlation based on glioma tissue analysis')+
  xlab("Order")+
  ylab("Correlation")
  

KUSTERcompare %>% dplyr::filter(colours == "High" & correlation > 0.5)
  



comparecor <- corKUSTER %>% 
  rename("correlation_Kuster" = "correlation") %>% 
  inner_join(corRNAprotNASPEAR %>% dplyr::select(gene, correlation), by = "gene") %>% 
  rename("correlation_us" = "correlation") %>% 
  dplyr::mutate(colours = ifelse(correlation_Kuster > 0.7 & correlation_us >0.7, "High", "Other")) %>% 
  dplyr::mutate(colours = ifelse((correlation_Kuster >0.3 & correlation_Kuster  <0.5) & (correlation_us >0.3 & correlation_us <0.5), "Moderate", colours))
  
pdf("corkuster.pdf", height = 4, width = 6)
ggplot(comparecor, aes( x= correlation_Kuster , y= correlation_us))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor(mapping = NULL)+
  labs (title = "Resemblance in correlation values for single genes", x= "Kuster's data" , y = "Our data")+
  theme_bw()

dev.off()
cor.test(comparecor$correlation_Kuster, comparecor$correlation_us)  


sum(comparecor$correlation_Kuster > 0.7 & comparecor$correlation_us >0.7)
sum((comparecor$correlation_Kuster >0.5 & comparecor$correlation_Kuster  <0.7) & (comparecor$correlation_us >0.5 & comparecor$correlation_us <0.7))
sum((comparecor$correlation_Kuster >0.3 & comparecor$correlation_Kuster  <0.5) & (comparecor$correlation_us >0.3 & comparecor$correlation_us <0.5))
sum((comparecor$correlation_Kuster >-0.3 & comparecor$correlation_Kuster  <0.3) & (comparecor$correlation_us >-0.3 & comparecor$correlation_us <0.3))

sum((comparecor$correlation_Kuster >-0.5 & comparecor$correlation_Kuster  < (-0.3)) & (comparecor$correlation_us >-0.5 & comparecor$correlation_us < (-0.3)))
sum(comparecor$correlation_Kuster < (-0.5) & comparecor$correlation_us < (-0.5))


#slangenplot met Kuster correlation
Slang_kuster_compare <- corrnatot %>% 
  dplyr::filter(Status == "Correlation") %>% 
  dplyr::select(gene, correlation, order1) %>% 
  rename("correlation_us" = "correlation") %>% 
  inner_join(corKUSTER %>% dplyr::select(gene, correlation), by = "gene") %>% 
  rename("correlation_Kuster" = "correlation") %>% 
  dplyr::mutate(colours = ifelse(correlation_Kuster > 0.7, "High Correlation", "0")) %>% 
  dplyr::mutate(colours = ifelse((correlation_Kuster >0.5 & correlation_Kuster  <0.7), "Moderate correlation", colours)) %>% 
  dplyr::mutate(colours = ifelse((correlation_Kuster >0.3 & correlation_Kuster  <0.5), "Low correlation", colours))


ggplot(Slang_kuster_compare, aes (x= order1, y= correlation_us, color = colours, label = gene)) +
  geom_point(data = subset(Slang_kuster_compare, colours == "0"), size = 2)+
  geom_point(data = subset(Slang_kuster_compare, colours == "Low correlation"), size = 1, shape = 4)+
  geom_point(data = subset(Slang_kuster_compare, colours == "Moderate correlation"), size = 1, shape = 4)+
  geom_point(data = subset(Slang_kuster_compare, colours == "High Correlation"), size = 1, shape = 4)+
  geom_text_repel(data = subset(Slang_kuster_compare, colours == "High Correlation", max.overlaps = 10))+
  theme_bw()+
  scale_color_manual(values=c("#008208","blue4", "yellow", "green")) +
  labs(title = "Correlation based on glioma tissue", color='Correlation among 27 tissues')+
  xlab("Order")+
  ylab("Correlation")
    