library(tidyr)
#Correlation with genelength
#Needs corRNAprotNASPEAR, gene.annot2, interprotRNA

tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))

gtf <- read.table('data/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf', comment.char='#',sep="\t",header=F)  %>% 
  dplyr::filter(V3 =='CDS') %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id ([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(protein_name = gsub("^.*protein_id([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(gene_name = gsub("^[ ]+","",gene_name))


gtf2 <- gtf %>% 
  dplyr::inner_join(interprotRNA)

#chr1	HAVANA	transcript	65419	71585	.	+	.	gene_id "ENSG00000186092.6"; transcript_id "ENST00000641515.2"; gene_type "protein_coding"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_name "OR4F5-202"; level 2; protein_id "ENSP00000493376.2"; hgnc_id "HGNC:14825"; tag "RNA_Seq_supported_partial"; tag "basic"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";

# load RNAseq ----
metadataRNA <- dplyr::left_join(gene.annot2, tmp, by = c("gene_id" = "Geneid"))%>% 
  dplyr::select(gene_id, gene_name, Chr, Start, Length, End, Strand)

metadataRNA <- dplyr::inner_join(interprotRNA, metadataRNA) %>% 
  dplyr::filter(gene_id %in% c("ENSG00000280987.4" , "ENSG00000187522.16","ENSG00000285053.1", "ENSG00000285441.1", "ENSG00000033050.9", "ENSG00000169093.16_PAR_Y", "ENSG00000169100.14_PAR_Y") == F)


lengthandcor <- left_join(corRNAprotNASPEAR, metadataRNA, by = c("gene" = "gene_name")) %>% 
  separate_rows(Chr:Strand, sep = ";")  %>% 
  dplyr::mutate(Exonlength = as.numeric(End)-as.numeric(Start))
  

grouped_lengthandcor<- lengthandcor %>%  group_by(gene)

maxexon <- grouped_lengthandcor %>% filter(Exonlength == max(as.numeric(Exonlength))) %>% unique()

maxlength <-grouped_lengthandcor %>% filter(End == max(as.numeric(End))) %>% distinct()

maxlength <-maxlength[!duplicated(maxlength$gene),]

minlength <- grouped_lengthandcor %>% filter(Start == min(as.numeric(Start)))

minlength <- minlength[!duplicated(minlength$gene),]

StartNew <- minlength %>%  dplyr::select(StartNew = Start)

chromosomelength <-  dplyr::left_join(maxlength, StartNew) %>% 
  dplyr::mutate(GeneLength = as.numeric(End)-as.numeric(StartNew) +1)

rm(StartNew)
                         
ggplot(chromosomelength, aes(x = log(1+GeneLength), y = correlation))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title = "Correlation increase with gene length", x= "log (1 + Gene length)" , y = "Correlation")+
  ggpubr::stat_cor(aes(shape = NULL, col= NULL, fill = NULL))+
  stat_poly_line()

pdf("corincreasegenelength.pdf", height = 4, width = 6)
ggplot(chromosomelength, aes(x = log(1+GeneLength), y = correlation))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_bw()+
  labs (title = "Correlation increase with gene length", x= "log (1 + Gene length)" , y = "Correlation")+
  ggpubr::stat_cor(aes(shape = NULL, col= NULL, fill = NULL))+
  stat_poly_line()

dev.off()

rm(StartNew, minlength, maxlength, maxexon, grouped_lengthandcor, gtf, tmp, gtf2, metadataRNA)
  


#slope(lengthandcor$Length, lengthandcor$correlation)


slope <- function(x,y) {
  return(cov(x,y)/var(x))
}


cor.test(x = log(1+chromosomelength$GeneLength), y = chromosomelength$correlation)



