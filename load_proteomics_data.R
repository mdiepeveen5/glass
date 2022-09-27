library("DEP")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("corrplot")

#load newdata met GS_ID names ----
#tsv tabel inladen
newdata <- read.table('~/projects/glass/data/glass/Proteomics/2022-03-31_data_update/20210729_084405_GLASSNL_LGG_dia_ProteinReport.tsv',
                      sep="\t"  ,quote="", header=T)

annotation_newdata <- readxl::read_xlsx('~/projects/glass//data//glass/Proteomics/2022-03-31_data_update/Annotation_Reduced_withControls.xlsx')

doublesProtein <- newdata %>% 
  dplyr::filter(PG.ProteinAccessions %in% c("P01892", "P04439", "P18463;P30475;Q29718;Q29836;Q95365","Q04826", "P01911", "Q29974","P42166","P42167")) 

#find the protein counts to exclude (ones with highest sumrow)
mostcountsdoubles <- doublesProtein %>% 
  dplyr::mutate_all(function(x){return(ifelse(x == "Filtered",1,0))}) %>% 
  dplyr::mutate(sumrow = rowSums(.))

#exclude 1,4,5,8
#exclude controls
newdata_doublesfiltered <- newdata %>% 
  dplyr::filter(PG.ProteinAccessions %in% c("P01892", "Q04826","P01911", "P42167") == F) %>% 
  dplyr::filter(PG.Genes != "")

rm(doublesProtein, mostcountsdoubles,newdata)
#newdata[1:6,8:10] %>%  dplyr::mutate_all(function(x){return(ifelse(x == "Filtered",1,0))})

#change file name to GLASS_ID
newdata_names <- newdata_doublesfiltered %>% 
  dplyr::mutate(PG.ProteinAccessions = NULL) %>% 
  dplyr::mutate(PG.CellularComponent = NULL) %>% 
  dplyr::mutate(PG.BiologicalProcess = NULL) %>% 
  dplyr::mutate(PG.MolecularFunction = NULL) %>%
  tibble::column_to_rownames("PG.Genes") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("File_Name_Proteomics") %>% 
  dplyr::mutate(File_Name_Proteomics = gsub('^.','', File_Name_Proteomics))

#NEW protein data on with GLASS ID + separated protein names
GLASS_proteomics <- annotation_newdata %>% 
  dplyr::select(File_Name_Proteomics, GS_ID) %>% 
  dplyr::left_join(newdata_names, by = "File_Name_Proteomics") %>% 
  dplyr::mutate(File_Name_Proteomics = NULL) %>% 
  dplyr::filter(GS_ID != "NA") %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("protein_name") %>% 
  separate_rows(protein_name, sep = ";")  

rm(annotation_newdata, newdata_doublesfiltered, newdata_names)

#load first protein and RNA data----

source("~/projects/glass/load_glass_data.R")

tmpprot <- read.csv('data/glass/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) %>% 
  dplyr::rename(Protein = X) %>% 
  dplyr::filter(Protein %in% c("HLA-B.1") == F)

tmpprot$Protein[tmpprot$Protein == "WARS"] <-"WARS1"
tmpprot$Protein[tmpprot$Protein == "NARS"] <-"NARS1"
tmpprot$Protein[tmpprot$Protein == "HARS"] <-"HARS1"
tmpprot$Protein[tmpprot$Protein == "DARS"] <-"DARS1"
tmpprot$Protein[tmpprot$Protein == "TARS"] <-"TARS1"
tmpprot$Protein[tmpprot$Protein == "VARS"] <-"VARS1"
tmpprot$Protein[tmpprot$Protein == "GARS"] <-"GARS1"
tmpprot$Protein[tmpprot$Protein == "IARS"] <-"IARS1"
tmpprot$Protein[tmpprot$Protein == "QARS"] <-"QARS1"
tmpprot$Protein[tmpprot$Protein == "AARS"] <-"AARS1"
tmpprot$Protein[tmpprot$Protein == "CARS"] <-"CARS1"
tmpprot$Protein[tmpprot$Protein == "SARS"] <-"SARS1"
tmpprot$Protein[tmpprot$Protein == "RARS"] <-"RARS1"
tmpprot$Protein[tmpprot$Protein == "YARS"] <-"YARS1"
tmpprot$Protein[tmpprot$Protein == "MARS"] <-"MARS1"
tmpprot$Protein[tmpprot$Protein == "KARS"] <-"KARS1"
tmpprot$Protein[tmpprot$Protein == "LARS"] <-"LARS1"

tmpprot$Protein[tmpprot$Protein == "HIST1H2BK"] <-"H2BC12"
tmpprot$Protein[tmpprot$Protein == "HIST1H2AB"] <-"H2AC4"
tmpprot$Protein[tmpprot$Protein == "HIST1H2AG"] <-"H2AC11"
tmpprot$Protein[tmpprot$Protein == "HIST1H1E"] <-"H1-4"
tmpprot$Protein[tmpprot$Protein == "HIST1H1B"] <-"H1-5"
tmpprot$Protein[tmpprot$Protein == "HIST1H1D"] <-"H1-3"
tmpprot$Protein[tmpprot$Protein == "HIST1H1C"] <-"H1-2"
tmpprot$Protein[tmpprot$Protein == "HIST1H2BO"] <-"H2BC17"
tmpprot$Protein[tmpprot$Protein == "HIST1H4A"] <-"H4C1"
tmpprot$Protein[tmpprot$Protein == "H3F3A"] <-"H3-3A"
tmpprot$Protein[tmpprot$Protein == "HIST2H2AC"] <-"H2AC20"
tmpprot$Protein[tmpprot$Protein == "HIST2H2AB"] <-"H2AC21"
tmpprot$Protein[tmpprot$Protein == "H1FX"] <-"H1-10"
tmpprot$Protein[tmpprot$Protein == "H2AFY2"] <-"MACROH2A2"
tmpprot$Protein[tmpprot$Protein == "H1F0"] <-"H1-0"
tmpprot$Protein[tmpprot$Protein == "H2AFY"] <-"MACROH2A1"
tmpprot$Protein[tmpprot$Protein == "H2AFZ"] <-"H2AZ1"

tmpprot$Protein[tmpprot$Protein == "TARSL2"] <-"TARS3"
tmpprot$Protein[tmpprot$Protein == "ASNA1"] <-"GET3"
tmpprot$Protein[tmpprot$Protein == "EPRS"] <-"EPRS1"
tmpprot$Protein[tmpprot$Protein == "ADSS"] <-"ADSS2"
tmpprot$Protein[tmpprot$Protein == "CRAD"] <-"CRACD"
tmpprot$Protein[tmpprot$Protein == "HIST1H2BA"] <-"H2BC1"
tmpprot$Protein[tmpprot$Protein == "KIF1BP"] <-"KIFBP"
tmpprot$Protein[tmpprot$Protein == "FAM129B"] <-"NIBAN2"
tmpprot$Protein[tmpprot$Protein == "FAM49A"] <-"CYRIA"
tmpprot$Protein[tmpprot$Protein == "RYDEN"] <-"SHFL"
tmpprot$Protein[tmpprot$Protein == "FAM49B"] <-"CYRIB"
tmpprot$Protein[tmpprot$Protein == "ADPRHL2"] <-"ADPRS"

#metadata tabel maken met patient en recurrent data en terugmatchen op RNA
metadata_prot_RNA <- data.frame(Sample_Name = colnames(tmpprot[,-1])) %>% 
  dplyr::mutate(sid = gsub('^[^_]+_[^_]+_([0-9]+)_.+$', '\\1', Sample_Name)) %>% 
  dplyr::mutate(resec = gsub('^[^_]+_[^_]+_[0-9]+_(.+)$', '\\1', Sample_Name)) %>% 
  dplyr::mutate(resec2 = case_when(
    resec == "P" ~ "R1",
    resec == "R1" ~ "R2",
    resec == "R2" ~ "R3",
    TRUE ~ "#$@#$"
  )) %>% 
  dplyr::mutate(Customer_ID=paste0(sid, "_", resec2)) %>% 
  dplyr::mutate(resec = NULL) %>% 
  dplyr::mutate(sid = NULL) %>%  
  dplyr::inner_join(metadata, by = "Customer_ID")

rm(metadata)

#Protein namen die verkeerd staan veranderen
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "WARS"] <-"WARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "NARS"] <-"NARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HARS"] <-"HARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "DARS"] <-"DARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "TARS"] <-"TARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "VARS"] <-"VARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "GARS"] <-"GARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "IARS"] <-"IARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "QARS"] <-"QARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "AARS"] <-"AARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "CARS"] <-"CARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "SARS"] <-"SARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "RARS"] <-"RARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "YARS"] <-"YARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "MARS"] <-"MARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "KARS"] <-"KARS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "LARS"] <-"LARS1"

GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H2BK"] <-"H2BC12"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H2AB"] <-"H2AC4"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H2AG"] <-"H2AC11"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H1E"] <-"H1-4"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H1B"] <-"H1-5"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H1D"] <-"H1-3"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H1C"] <-"H1-2"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H2BO"] <-"H2BC17"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H4A"] <-"H4C1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H3F3A"] <-"H3-3A"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST2H2AC"] <-"H2AC20"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST2H2AB"] <-"H2AC21"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H1FX"] <-"H1-10"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H2AFY2"] <-"MACROH2A2"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H1F0"] <-"H1-0"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H2AFY"] <-"MACROH2A1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "H2AFZ"] <-"H2AZ1"

GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "TARSL2"] <-"TARS3"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "ASNA1"] <-"GET3"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "EPRS"] <-"EPRS1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "ADSS"] <-"ADSS2"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "CRAD"] <-"CRACD"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "HIST1H2BA"] <-"H2BC1"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "KIF1BP"] <-"KIFBP"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "FAM129B"] <-"NIBAN2"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "FAM49A"] <-"CYRIA"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "RYDEN"] <-"SHFL"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "FAM49B"] <-"CYRIB"
GLASS_proteomics$protein_name[GLASS_proteomics$protein_name == "ADPRHL2"] <-"ADPRS"

#gene annotations for RNA data----
gene.annot <- read.table('data/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf', comment.char='#',sep="\t",header=F)
gene.annot2 <- gene.annot %>% 
  dplyr::filter(V3 =='gene') %>% 
  dplyr::mutate(gene_id = gsub("^.*gene_id ([^;]+);.*$","\\1",V9)) %>% 
  dplyr::mutate(gene_name = gsub("^.*gene_name([^;]+);.*$","\\1",V9)) %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::mutate(gene_name = gsub("^[ ]+","",gene_name))

rm(gene.annot)

#RNA filteren op patienten en goed benamen----
#patienten die in RNA en protein zitten
RNA_counts_filtered <- readcount.vst %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("patient_id") %>% 
  dplyr::inner_join(metadata_prot_RNA %>% dplyr::select(GS_ID), by = c("patient_id" ="GS_ID")) %>% 
  tibble::column_to_rownames("patient_id") %>% 
  t() 

readcounts_filtered <- readcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("patient_id") %>% 
  dplyr::inner_join(metadata_prot_RNA %>% dplyr::select(GS_ID), by = c("patient_id" ="GS_ID")) %>% 
  tibble::column_to_rownames("patient_id") %>% 
  t() 

rm(readcount.vst, readcounts)

#get out double read of RNA
#doublesRNA <- RNA_counts_filtered[c("ENSG00000280987.4", "ENSG00000015479.19", "ENSG00000284024.2", "ENSG00000187522.16", "ENSG00000285053.1", "ENSG00000284770.2"),] %>% 
#  as.data.frame() %>% 
#  dplyr::mutate(sumrow = rowSums(.))

#filter dubbele genen in RNA met laagste counts eruit
RNA_counts_filtered <-RNA_counts_filtered %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::filter(gene_id %in% c("ENSG00000280987.4" , "ENSG00000187522.16","ENSG00000285053.1") == F)

readcounts_filtered <-readcounts_filtered %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('gene_id') %>% 
  dplyr::filter(gene_id %in% c("ENSG00000280987.4" , "ENSG00000187522.16","ENSG00000285053.1") == F)

#genen in protein gefilterd uit gene.annot2
corresponding_genes <- data.frame(tmpprot$Protein) %>% 
  dplyr::inner_join(gene.annot2, by = c("tmpprot.Protein" = "gene_name")) %>% 
  as.data.frame()

#RNA with proper gene names
total_RNA <- corresponding_genes %>% 
  dplyr::inner_join(RNA_counts_filtered) %>% 
  dplyr::mutate(gene_id =NULL) %>% 
  dplyr::rename(gene_name = tmpprot.Protein)

readcounts_filteredtotal <- corresponding_genes %>% 
  dplyr::inner_join(readcounts_filtered) %>% 
  dplyr::mutate(gene_id =NULL) %>% 
  dplyr::rename(gene_name = tmpprot.Protein)

rm(RNA_counts_filtered, readcounts_filtered)

#protein data filteren op patienten in RNA en protein----
pre_protein_counts_filtered <- tmpprot %>% 
  tibble::column_to_rownames("Protein") %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("patient_id")

#GS_ID veranderen en filteren op genen
protein_counts_filtered <- dplyr::select(metadata_prot_RNA, "GS_ID", "Sample_Name.x") %>% 
  left_join(pre_protein_counts_filtered, by = c("Sample_Name.x" = "patient_id")) %>% 
  dplyr::mutate(Sample_Name.x=NULL) %>% 
  tibble::column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Protein") %>% 
  dplyr::inner_join(corresponding_genes , by = c("Protein" = "tmpprot.Protein" )) %>% 
  dplyr::mutate(gene_id=NULL) %>% 
  distinct() %>% 
  as.data.frame()

rm(pre_protein_counts_filtered)

#filteren op genen in RNA en protein
interprotRNA <- data.frame(gene_name = intersect(protein_counts_filtered$Protein, total_RNA$gene_name))


protein_counts_filtered <- dplyr::inner_join(protein_counts_filtered, interprotRNA, by = c("Protein" = "gene_name"))
total_RNA <- dplyr::left_join(interprotRNA, total_RNA, by = "gene_name")

rm(tmpprot)


#filter nieuwe protein data op proteins die in RNA en protein zitten + separate protein names
GLASS_interprotRNA <- GLASS_proteomics%>% 
  right_join(interprotRNA, by = c("protein_name" ="gene_name"))


rm(GLASS_proteomics)

#lijst met patienten die in beide datasets zitten
GS_ID_RNA_protein <- total_RNA %>% 
  tibble::column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GS_ID") %>% 
  dplyr::select(GS_ID)

#raw protein data filtered on GS_ID in RNA and protein
GLASS_filtered <- GLASS_interprotRNA %>% 
  column_to_rownames("protein_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GS_ID") %>% 
  dplyr::right_join(GS_ID_RNA_protein) %>% 
  column_to_rownames("GS_ID") %>% 
  t() %>% 
  as.data.frame()

rm(GLASS_interprotRNA)

filteredNArawprotein <- as.data.frame(GLASS_filtered, stringsAsFactors = FALSE) %>% 
  dplyr::mutate_all(function(x){return(ifelse(x == "Filtered",NA_integer_,x))})

filteredNArawprotein[] <- lapply(filteredNArawprotein, function (x) as.numeric(as.character(x)))




#check how many filtered entries a gene has:
#filteredentries_rawprotein <- GLASS_filtered %>% 
#  dplyr::mutate_all(function(x){return(ifelse(x == "Filtered",1,0))}) %>% 
#  dplyr::mutate(sumrow = rowSums(.)) %>% 
#  dplyr::filter(sumrow<4) %>% 
#  rownames_to_column("protein_name") %>% 
#  dplyr::select(protein_name)