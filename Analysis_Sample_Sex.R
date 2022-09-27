source("~/projects/glass/load_glass_data.R")

# supervised analysis ----

## deseq2 ----

#stopifnot(colnames(totaldata) == coldata$GS_ID)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts,
                                      colData = metadata,
                                      design= ~Sample_Sex)
dds <-DESeq(dds)
res <-results(dds) %>% 
  as.data.frame %>% 
  dplyr::arrange(padj,pvalue)

selected.features.supervised <- res %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::select('gene_id') %>% 
  dplyr::slice_head(n=1000)

readcount.vst <-readcount.vst %>% 
  tibble::rownames_to_column('gene_id')


supervised <- dplyr::left_join(selected.features.supervised, readcount.vst, by = "gene_id") %>% 
  tibble::column_to_rownames('gene_id')

pca_supervised <- supervised %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')

supervised_sorted <-dplyr::left_join(metadata, pca_supervised)

ggplot(supervised_sorted, aes(x=PC1, y=PC2, color = Sample_Sex, group = GLASS_ID)) +
  geom_point() +
  geom_line() +
  theme_classic()

readcount.vst <-readcount.vst %>% 
  tibble::column_to_rownames('gene_id')

##volcanoplot ----

gene.annot <- read.table('data/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf', comment.char='#',sep="\t",header=F)
gene.annot <- gene.annot %>% 
  dplyr::filter(V3 =='gene') %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id ([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name([^;]+);.*$","\\1",V9)) %>% 
  dplyr::select(gene_id, gene_name)


gene_annotations <- res %>%
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::left_join(gene.annot, by = 'gene_id') %>% 
  tibble::column_to_rownames('gene_id')


library(EnhancedVolcano)
EnhancedVolcano(gene_annotations[1:6],
                lab = gene_annotations$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue')

#highRNAexpres <- left_join(select(gene_annotations$gene_name), selected.features.supervised$gene_id)

# unsupervised analysis ----
selected.features <- readcount.vst %>% 
  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>% 
  dplyr::arrange(-mad) %>% 
  dplyr::slice_head(n=1000) %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::select('gene_id')

readcount.vst <-readcount.vst %>% 
  tibble::rownames_to_column('gene_id')
  
total <- dplyr::left_join(selected.features, readcount.vst, by = "gene_id") %>% 
  tibble::column_to_rownames('gene_id')
  

#gen namen sorteden op padj
pca <- total %>% 
  t %>% 
  prcomp() %>% 
  purrr::pluck('x') %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column('GS_ID')

unsupervised_sorted <-dplyr::left_join(metadata, pca)

ggplot(unsupervised_sorted, aes(x=PC3, y=PC4, color = Sample_Sex)) + geom_point()
