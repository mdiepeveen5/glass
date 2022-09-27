#needs; interprotRNA

tmp <- read.delim('data/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T)
colnames(tmp) <- gsub(".Aligned.+bam$","",gsub("^X.+new.","",colnames(tmp)))

unfiltered_RNA <- dplyr::left_join(gene.annot2, tmp, by = c("gene_id" = "Geneid")) %>% 
  dplyr::filter(gene_id %in% c("ENSG00000169100.14_PAR_Y", "ENSG00000280987.4", "ENSG00000285053.1","ENSG00000187522.16", "ENSG00000169093.16_PAR_Y", "ENSG00000285441.1","ENSG00000033050.9")==F,)


unfiltered_protein <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) %>% 
  dplyr::rename(Protein = X)

pre_protein_counts_filtered <- unfiltered_protein %>% 
  tibble::column_to_rownames("Protein") %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("patient_id")

protein_counts_filtered_allprotein <- dplyr::select(metadata_prot_RNA, "GS_ID", "Sample_Name.x") %>% 
  left_join(pre_protein_counts_filtered, by = c("Sample_Name.x" = "patient_id")) %>% 
  dplyr::mutate(Sample_Name.x=NULL) %>% 
  tibble::column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Protein")

rm(pre_protein_counts_filtered)

protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "WARS"] <-"WARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "NARS"] <-"NARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HARS"] <-"HARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "DARS"] <-"DARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "TARS"] <-"TARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "VARS"] <-"VARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "GARS"] <-"GARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "IARS"] <-"IARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "QARS"] <-"QARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "AARS"] <-"AARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "CARS"] <-"CARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "SARS"] <-"SARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "RARS"] <-"RARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "YARS"] <-"YARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "MARS"] <-"MARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "KARS"] <-"KARS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "LARS"] <-"LARS1"

protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H2BK"] <-"H2BC12"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H2AB"] <-"H2AC4"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H2AG"] <-"H2AC11"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H1E"] <-"H1-4"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H1B"] <-"H1-5"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H1D"] <-"H1-3"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H1C"] <-"H1-2"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H2BO"] <-"H2BC17"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H4A"] <-"H4C1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H3F3A"] <-"H3-3A"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST2H2AC"] <-"H2AC20"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST2H2AB"] <-"H2AC21"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H1FX"] <-"H1-10"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H2AFY2"] <-"MACROH2A2"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H1F0"] <-"H1-0"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H2AFY"] <-"MACROH2A1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "H2AFZ"] <-"H2AZ1"



protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "TARSL2"] <-"TARS3"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "ASNA1"] <-"GET3"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "EPRS"] <-"EPRS1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "ADSS"] <-"ADSS2"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "CRAD"] <-"CRACD"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "HIST1H2BA"] <-"H2BC1"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "KIF1BP"] <-"KIFBP"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "FAM129B"] <-"NIBAN2"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "FAM49A"] <-"CYRIA"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "RYDEN"] <-"SHFL"
protein_counts_filtered_allprotein$Protein[protein_counts_filtered_allprotein$Protein == "FAM49B"] <-"CYRIB"



unfiltered_RNA2<- inner_join(protein_counts_filtered_allprotein, unfiltered_RNA, by = c("Protein" = "gene_name"))

protnotRNA<- anti_join(protein_counts_filtered_allprotein, unfiltered_RNA2, by = c("Protein"= "Protein")) %>% 
  dplyr::select(Protein)


rm(tmp, unfiltered_RNA, all_data_unfiltered)


print(as.data.frame(protnotRNA), row.names = FALSE)


print(as.data.frame(protein_counts_filtered_allprotein$Protein[2994:3250],), row.names = F)
