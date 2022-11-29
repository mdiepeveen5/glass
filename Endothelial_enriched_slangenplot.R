#Slangenplot endothelial

library(readr)
endothelial_cells <- read_csv("endothelial_cells.csv", show_col_types = FALSE) %>% 
  dplyr::select(name)

slang_endothelial <- corrnatot %>%  dplyr::filter(Status == "Correlation") %>% 
  dplyr::select(gene, correlation, order1) %>% 
  dplyr::mutate(colors = ifelse(gene %in% endothelial_cells$name, "Endothelial cell", "Other"))

meanendo <- slang_endothelial %>% dplyr::filter(colors == "Endothelial cell") %>% dplyr::select(correlation)
mean(meanendo$correlation)

meanother <- slang_endothelial %>% dplyr::filter(colors == "Other") %>% dplyr::select(correlation)
mean(meanother$correlation)

t.test(meanendo$correlation, meanother$correlation, var.equal = T)



ggplot(slang_endothelial, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(slang_endothelial, colors == "Other"), size = 1)+
  geom_point(data = subset(slang_endothelial, colors == "Endothelial cell"), size = 1)+
  geom_text_repel(data = subset(slang_endothelial, colors == "Endothelial cell", max.overlaps = 10))+
  scale_color_manual(values=c("red","blue4"))+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")


rm(slang_endothelial, meanendo, meanother)


molecular_layer <- read_csv("molecular_layer.csv", show_col_types = FALSE) %>% 
  dplyr::select(name)


slang_molecular_layer <- corrnatot %>%  dplyr::filter(Status == "Correlation") %>% 
  dplyr::select(gene, correlation, order1) %>% 
  dplyr::mutate(colors = ifelse(gene %in% molecular_layer$name, "Molecular layer cell", "Other"))

meanmol <- slang_molecular_layer %>% dplyr::filter(colors == "Molecular layer cell") %>% dplyr::select(correlation)
mean(meanmol$correlation)

meanothermol <- slang_molecular_layer %>% dplyr::filter(colors == "Other") %>% dplyr::select(correlation)
mean(meanothermol$correlation)

t.test(meanmol$correlation, meanothermol$correlation,  var.equal = T)



ggplot(slang_molecular_layer, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(slang_molecular_layer, colors == "Other"), size = 1)+
  geom_point(data = subset(slang_molecular_layer, colors == "Molecular layer cell"), size = 1)+
  geom_text_repel(data = subset(slang_molecular_layer, colors == "Molecular layer cell", max.overlaps = 10))+
  scale_color_manual(values=c("red","blue4"))+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")

rm(slang_molecular_layer, meanendo, meanother)



neuropil <- read_csv("neuropil.csv", show_col_types = FALSE) %>% 
  dplyr::select(name)

slang_neuropil <- corrnatot %>%  dplyr::filter(Status == "Correlation") %>% 
  dplyr::select(gene, correlation, order1) %>% 
  dplyr::mutate(colors = ifelse(gene %in% neuropil$name, "Neuropil", "Other"))


meanneuropil <- slang_neuropil %>% dplyr::filter(colors == "Neuropil") %>% dplyr::select(correlation)
mean(meanmol$correlation)

meanotherneuropil <- slang_neuropil %>% dplyr::filter(colors == "Other") %>% dplyr::select(correlation)
mean(meanothermol$correlation)


t.test(meanneuropil$correlation, meanotherneuropil$correlation, var.equal = T)



ggplot(slang_neuropil, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(slang_neuropil, colors == "Other"), size = 1)+
  geom_point(data = subset(slang_neuropil, colors == "Neuropil"), size = 1)+
  geom_text_repel(data = subset(slang_neuropil, colors == "Neuropil", max.overlaps = 10))+
  scale_color_manual(values=c("red","blue4"))+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")

rm(slang_neuropil, meanendo, meanother)



all_enriched <- slang_endothelial %>% 
  dplyr::mutate(colors = ifelse(gene %in% molecular_layer$name, "Molecular layer cell", colors)) %>% 
  dplyr::mutate(colors = ifelse(gene %in% neuropil$name, "Neuropil", colors)) %>% 
  dplyr::mutate(colors = ifelse(gene %in% granular_cells$name, "Granular Cell", colors))  
  

  
ggplot(all_enriched, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(all_enriched, colors == "Other"), size = 1)+
  geom_point(data = subset(all_enriched, colors == "Endothelial cell"), size = 1)+
  geom_point(data = subset(all_enriched, colors == "Molecular layer cell"), size = 1)+
  geom_point(data = subset(all_enriched, colors == "Neuropil"), size = 1)+
  geom_point(data = subset(all_enriched, colors == "Granular Cell"), size = 1)+
  scale_color_manual(values=c("yellow", "red","blue4", "green", "grey"))+
  theme_bw()+
  labs(title = "Correlation", color='Cell type')+
  xlab("Order")+
  ylab("Correlation")

rm(all_enriched, endothelial_cells, molecular_layer, neuropil, slang_endothelial, slang_molecular_layer, slang_neuropil, meanother, meanothermol, meanneuropil, meanothergranular, meanotherneuropil, meanendo, meangranular, meanmol, slang_granular)


granular_cells <- read_csv("granular_cells.csv")

slang_granular <- corrnatot %>%  dplyr::filter(Status == "Correlation") %>% 
  dplyr::select(gene, correlation, order1) %>% 
  dplyr::mutate(colors = ifelse(gene %in% granular_cells$name, "Granular Cell", "Other"))


meangranular <- slang_granular %>% dplyr::filter(colors == "Granular Cell") %>% dplyr::select(correlation)
mean(meangranular$correlation)

meanothergranular <- slang_granular %>% dplyr::filter(colors == "Other") %>% dplyr::select(correlation)
mean(meanothergranular$correlation)


t.test(meangranular$correlation, meanothergranular$correlation, var.equal = T)



ggplot(slang_granular, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(slang_granular, colors == "Other"), size = 1)+
  geom_point(data = subset(slang_granular, colors == "Neuropil"), size = 1)+
  geom_text_repel(data = subset(slang_granular, colors == "Neuropil", max.overlaps = 10))+
  scale_color_manual(values=c("red","blue4"))+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")


##MEAN VALUES PLOT



 
slang_neuropil$colors <- factor(slang_neuropil$colors, levels = c("Neuropil", "Other"))
 
pdf("slang_neuropil.pdf", height = 4, width = 4)

ggboxplot(slang_neuropil, x = "colors", y = "correlation",
          color = "colors", palette = c("#00AFBB", "#E7B800"))+
  stat_compare_means(method = "t.test", paired = F)

dev.off()
t.test(meanneuropil$correlation, meanotherneuropil$correlation)

compare_means(correlation~colors, data = slang_neuropil, method = "t.test")

slang_granular$colors <- factor(slang_granular$colors, levels = c("Granular Cell", "Other"))

pdf("slang_granular.pdf", height = 4, width = 4)

ggboxplot(slang_granular, x = "colors", y = "correlation",
          color = "colors", palette = c("#00AFBB", "#E7B800"))+
  stat_compare_means(method = "t.test", paired = F)
dev.off()

compare_means(correlation~colors, data = slang_granular, method = "t.test")


slang_molecular_layer$colors <- factor(slang_molecular_layer$colors, levels = c("Molecular layer cell", "Other"))

pdf("slang_molecular_layer.pdf", height = 4, width = 4)

ggboxplot(slang_molecular_layer, x = "colors", y = "correlation",
          color = "colors", palette = c("#00AFBB", "#E7B800"))+
  stat_compare_means(method = "t.test", paired = F)
dev.off()


compare_means(correlation~colors, data = slang_molecular_layer, method = "t.test")


slang_endothelial$colors <- factor(slang_endothelial$colors, levels = c("Endothelial cell", "Other"))

pdf("slang_endothelial.pdf", height = 4, width = 4)
ggboxplot(slang_endothelial, x = "colors", y = "correlation",
          color = "colors", palette = c("#00AFBB", "#E7B800"))+
  stat_compare_means(method = "t.test", paired = F)


dev.off()
compare_means(correlation~colors, data = slang_molecular_layer, method = "t.test")

  