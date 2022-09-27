##Figures for paper
load(file = "corrnatotPAPER.Rdata")

#with non imputed values



corrnatotPAPER$group <- "0"
corrnatotPAPER$group[corrnatotPAPER$gene %in% c("LMNB1", "ACTC1", "CDK1", "H4C1", "H2AC4", "H2AC11", "H2BC12", "H1-5", "H2BC17")] <- "up.1"
corrnatotPAPER$group[corrnatotPAPER$gene %in% c("COL6A3", "COL6A2", "TGFBI","HSPG2","LUM","FCGBP","CHI3L1","COL1A2","ACAN","VGF","COL4A2","COL1A1")] <- "up.2"
#corrnatot$group[corrnatot$gene %in% c("AQP1", "FABP5")] <- "up.3"
corrnatotPAPER$group[corrnatotPAPER$gene %in% c("MPP6", "OGFRL1", "GJA1", "RGS6")] <- "down"

corrnatotPAPER <- corrnatotPAPER %>%  dplyr::mutate(colors = ifelse(group == "0", Status, group)) %>% 
  dplyr::mutate(colors = ifelse(Status == "Random Correlation" & colors != "Random Correlation", Status, colors))

ggplot(corrnatotPAPER, aes (x= order1, y= correlation, color = colors, label = gene)) +
  geom_point(data = subset(corrnatotPAPER, group == '0'), size=2) +
  geom_point(data = subset(corrnatotPAPER, group != '0'), size=2) +
  geom_text_repel(data = subset(corrnatotPAPER, group != '0' & Status == "Correlation" )) +
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")


#inkscape

ggplot(corrnatotPAPER, aes (x= order1, y= correlation, color = Status))+
  geom_point()+
  theme_bw()+
  labs(title = "Correlation vs Random Correlation", color='Correlation and cluster type')+
  xlab("Order")+
  ylab("Correlation")

ggsave(file = "coorrnatotPAPERplot.pdf", plot = corrnatotPAPERplot)


mean(corrnatotPAPER$correlation %>% dplyr::filter(Status == "Correlation"))

corrna1 <- corrnatotPAPER %>%  dplyr::filter(Status == "Correlation") 
mean(corrna1$correlation)

corrna2 <- corrnatotPAPER %>%  dplyr::filter(Status == "Random Correlation")
mean(corrna2$correlation)

wilcox.test(corrna1$correlation, corrna2$correlation)


last100 <- corrna1 %>% dplyr::filter(order1 > 3083) 
print(as.data.frame(last100$gene), row.names = F)

first100 <- corrna1 %>% dplyr::filter(order1 <101)
