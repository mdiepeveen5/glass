differential_RNA <- deseq2Results %>% 
  left_join(gene.annot2)

signupRNA <- differential_RNA %>% 
  dplyr::filter(log2FoldChange >0.75 & padj <0.05)

signdownRNA <- differential_RNA %>% 
  dplyr::filter(log2FoldChange < (-0.75) & padj <0.05)




write.table(as.data.frame(signupRNA$gene_name), file = "txtsignup.txt", sep = "\t", row.names = F)


print(as.data.frame(signupRNA$gene_name), row.names = FALSE)
print(as.data.frame(signdownRNA$gene_name), row.names = FALSE)
