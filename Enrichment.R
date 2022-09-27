library(Hmisc)

res.paired.a.exon <- readRDS("~/projects/glass/res.paired.a.exon.Rds")

# pathway enrichment
#https://biit.cs.ut.ee/gplink/l/CDRGUemNQi
#https://biit.cs.ut.ee/gplink/l/HTRVSdREQE    all proteins

# g:profiler
# co-expressie <rna/protein>
#fgsea
rownames(corRNAprothigh) <- NULL

print(as.data.frame(corRNAprothigh[1:197,]), row.names = FALSE)

corabove0.6 <- corRNAprotNASPEAR %>% 
  dplyr::filter(correlation >0.5)

#GSEA----
# ALL pathways
pathways.MSigDB <- fgsea::gmtPathways('msigdb.v7.5.1.entrez.gmt')
pathways.MSigDB <- data.frame(ENTREZ = unlist(pathways.MSigDB), SubCategory = rep(names(pathways.MSigDB), S4Vectors::elementNROWS(pathways.MSigDB)))

# Get ENSEMBL ids.
geneIDs <- clusterProfiler::bitr(pathways.MSigDB$ENTREZ, fromType = 'ENTREZID', toType = 'ENSEMBL', OrgDb = org.Hs.eg.db::org.Hs.eg.db)
geneIDs <- geneIDs[!duplicated(geneIDs$ENTREZID),]

#genes are filtered on highest genecount 
gseaGenesCOR <- tibble(corRNAprothigh) %>% dplyr::inner_join(tibble(corresponding_genes), by =c("gene" = "tmpprot.Protein")) %>% 
  dplyr::filter(gene_id %in% c("ENSG00000169100.14_PAR_Y", "ENSG00000280987.4", "ENSG00000285053.1","ENSG00000187522.16", "ENSG00000169093.16_PAR_Y", "ENSG00000285441.1","ENSG00000033050.9")==F,) %>% 
  dplyr::mutate(newgeneid = gsub("\\..*","", gene_id)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL")) %>% 
  dplyr::mutate(rank = 1:n()) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000185596", "374666", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000214338", "387104", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000277203", "8263", ENTREZID)) 

gseaGenesRNA <- res.paired.a.exon %>% 
  dplyr::mutate(newgeneid = gsub("_.*","", gene_uid)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL")) %>% 
  dplyr::mutate(rank = 1:n())
  

pws <- fgsea::gmtPathways('msigdb.v7.5.1.entrez.gmt')

pwshall <- fgsea::gmtPathways('h.all.v7.5.1.entrez.gmt')

#make the ranks
ranks <- gseaGenesCOR$rank
names(ranks) = gseaGenesCOR$ENTREZID 

ranksRNA <- gseaGenesRNA$rank
names(ranksRNA) = gseaGenesRNA$ENTREZID

#GSEA analysis
#pathways die bovenaan enrichment staan opzoeken 
fgseaResCOR <-fgsea(pathways = pws, stats = ranks)
fgseaResCOR <- fgseaResCOR[order(fgseaResCOR$padj),]

fgseaResCORhall <-fgsea(pathways = pwshall, stats = ranks)
fgseaResCORhall <- fgseaResCORhall[order(fgseaResCORhall$pval)]

fgseaResRNA <- fgsea(pathways = pws, stats = ranksRNA)
fgseaResRNA <- fgseaResRNA[order(fgseaResRNA$padj),]

fgseaResRNAhall <-fgsea(pathways = pwshall, stats = ranksRNA)
fgseaResRNAhall <-fgseaResRNAhall[order(fgseaResRNAhall$padj),]


plotEnrichment(pws$"KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN", stats = ranksRNA, gseaParam = 1, ticksSize = 0.2)

dev.off()

plotEnrichment(pwshall$"HALLMARK_HEME_METABOLISM", stats = ranks, gseaParam = 1, ticksSize = 0.2)

#plot GSEA table

toppwsCOR <- fgseaResCOR[head(order(pval), n = 15)][order(NES), pathway]

plotGseaTable(pws[toppwsCOR], ranks, fgseaResCOR, gseaParam = 0.5)

toppwsRNA <- fgseaResRNA[head(order(pval), n = 15)][order(NES), pathway]

plotGseaTable(pws[toppwsRNA], ranksRNA, fgseaResRNA, gseaParam = 0.5)


#generate pathway met clustergenen van Youri, 4 pathways, up1 up2 up3 down --> zelfde structuur maken als pws met entez id's 

#pathways van youri maken en entrezid's toevoegen
up1 = dge.partially.paired.clusters %>% 
  dplyr::filter(up.1 == T) %>% 
  dplyr::select(gene_name, gene_uid) %>% 
  dplyr::mutate(newgeneid = gsub("_.*","", gene_uid)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL")) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000258910", "107985816", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000237380", "100506783", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000249362", "111082994", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000253819", "104266958", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000285077", "89839", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000245694", "105378481", ENTREZID))  

up2 = dge.partially.paired.clusters %>% 
  dplyr::filter(up.2 == T) %>% 
  dplyr::select(gene_name, gene_uid)%>% 
  dplyr::mutate(newgeneid = gsub("_.*","", gene_uid)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL"))%>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000273143", "105378481", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000257139", "105369826", ENTREZID))

up3 = dge.partially.paired.clusters %>% 
  dplyr::filter(up.3 == T) %>% 
  dplyr::select(gene_name, gene_uid)%>% 
  dplyr::mutate(newgeneid = gsub("_.*","", gene_uid)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL")) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000237166", "100507140", ENTREZID)) %>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000177822", "90768", ENTREZID))%>% 
  dplyr::mutate(ENTREZID = ifelse(newgeneid == "ENSG00000259380", "120017343", ENTREZID))

down = dge.partially.paired.clusters %>% 
  dplyr::filter(down == T) %>% 
  dplyr::select(gene_name, gene_uid)%>% 
  dplyr::mutate(newgeneid = gsub("_.*","", gene_uid)) %>% 
  dplyr::left_join(geneIDs, by = c("newgeneid" = "ENSEMBL"))

pathup1 = list('pathup1' = up1$ENTREZID)
pathup2 = list ('pathup2' = up2$ENTREZID)
pathup3 = list ('pathup3' = up3$ENTREZID)
pathdown = list('pathdown' = down$ENTREZID)

pws_appended <-pathup1 %>% 
  append(pathup2) %>% 
  append(pathup3) %>% 
  append(pathdown)

#protein RNA intersection 
interup1 <- up1 %>% 
  dplyr::inner_join(interprotRNA)

interup2 <- up2 %>% 
  dplyr::inner_join(interprotRNA)

interup3 <- up3 %>% 
  dplyr::inner_join(interprotRNA)

interdown <- down %>% 
  dplyr::inner_join(interprotRNA)

interpathup1 = list('interpathup1' = interup1$ENTREZID)
interpathup2 = list ('interpathup2' = interup2$ENTREZID)
interpathup3 = list ('interpathup3' = interup3$ENTREZID)
interpathdown = list('interpathdown' = interdown$ENTREZID)


pws_inter <- interpathup1 %>% 
  append(interpathup2) %>% 
  append(interpathup3) %>% 
  append(interpathdown)

fgseaResaddedpws <-fgsea(pathways = pws_inter, stats = ranks)
fgseaResaddedpws <- fgseaResaddedpws[order(fgseaResaddedpws$padj)]


#lijst van youri met van alle genen p-value, daarmee pathway enrichment doen
#die vergelijken met pathway enrichment van correlatie

#g:profiler

wilcox.test(corRNAprotNA$correlation)

rm(pathup1, pathup2, pathup3, interpathdown, interpathup1, interpathup2, interpathup3)