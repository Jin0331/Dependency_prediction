library(R.utils)
library(biomaRt)
library(tidyverse)
library(GEOquery)
library(minfi)
library(RMariaDB)
library(progress)
library(Prep4DeepDEP)

# db connection
con <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!",
                      dbname = "CCLE_22Q2")  

# TCGA PANCAN for pretrain ####
## raw data load
{
  tcga_exp_raw <- data.table::fread(file = "RAW/PANCANCER/tcga_RSEM_gene_tpm.gz") %>% 
    as_tibble()
  
  variant_type <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site")
  tcga_mut_raw <- data.table::fread(file = "RAW/PANCANCER/mc3.v0.2.8.PUBLIC.maf.gz") %>% 
    as_tibble() %>% 
    filter(Variant_Classification %in% variant_type) %>% 
    separate(col = Tumor_Sample_Barcode, into = LETTERS[1:7], sep = "-") %>% 
    select(-E,-F,-G) %>% 
    unite(col = "Tumor_Sample_Barcode", A,B,C,D, sep = "-") %>% 
    filter(str_detect(string = Tumor_Sample_Barcode, pattern = "^*A$")) %>% 
    mutate(Tumor_Sample_Barcode =  substring(Tumor_Sample_Barcode, 1, nchar(Tumor_Sample_Barcode) - 1))
  
  tcga_cna_raw <- data.table::fread(file = "RAW/PANCANCER/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg") %>% 
    as_tibble() %>% 
    dplyr::rename(Tumor_Sample_Barcode = Sample) %>% 
    separate(col = Tumor_Sample_Barcode, into = LETTERS[1:7], sep = "-") %>% 
    select(-E,-F,-G) %>% 
    unite(col = "Tumor_Sample_Barcode", A,B,C,D, sep = "-") %>% 
    filter(str_detect(string = Tumor_Sample_Barcode, pattern = "^*A$")) %>% 
    mutate(Tumor_Sample_Barcode =  substring(Tumor_Sample_Barcode, 1, nchar(Tumor_Sample_Barcode) - 1))
  
  tcga_meth_raw <- data.table::fread(file = "RAW/PANCANCER/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv") %>% 
    as_tibble() %>% 
    dplyr::rename(Probe = `Composite Element REF`)
}

if(!file.exists("TCGA_OMICS_INTER_SAMPLE_BARCODE.txt")){
  tcga_primary_sample <- read_delim(file = "RAW/PANCANCER/TCGA_phenotype_denseDataOnlyDownload.tsv",
                                    delim = "\t") %>% 
    filter(sample_type == "Primary Tumor") %>% 
    pull(sample)
  
  ## exp
  tcga_exp_sample <- tcga_exp_raw %>% 
    as_tibble() %>% 
    colnames() %>% unique()
  
  ## mut
  tcga_mut_sample <- tcga_mut_raw %>% 
    pull(Tumor_Sample_Barcode) %>% 
    unique()
  
  # cna
  tcga_cna_sample <- tcga_cna_raw %>% 
    pull(Tumor_Sample_Barcode) %>% 
    unique() 
  
  # methlation
  tcga_meth_sample <- tcga_meth_raw %>% colnames() %>% .[-1] %>% 
    lapply(X = ., FUN = function(s){
      tmp <- s %>% str_split(pattern = "-") %>% unlist() %>% 
        .[1:4] %>% 
        paste0(collapse = "-") %>% 
        substring(., 1, nchar(.) - 1)
      return(tmp)
    }) %>% unlist() %>% unique()
  colnames(tcga_meth_raw) <- c("Probe", tcga_meth_sample)
  
  tcga_omics_sample <- intersect.Vector(tcga_exp_sample, tcga_mut_sample) %>% 
    intersect.Vector(., tcga_cna_sample) %>% 
    intersect.Vector(., tcga_meth_sample)
  
  # p_tcga <- ggVennDiagram::ggVennDiagram(x = list(Exp = tcga_exp_sample,
  #                                                 Mutation = tcga_mut_sample,
  #                                                 Cna = tcga_cna_sample,
  #                                                 Meth = tcga_meth_sample)) +
  #   ggtitle("TCGA-Omics integration") + 
  #   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  #   theme(legend.position = "none")
  
  tcga_omics_sample %>% tibble(Tumor_Sample_Barcode = .) %>% 
    write_delim(file = "TCGA_OMICS_INTER_SAMPLE_BARCODE.txt", delim = '\t')
  ggsave(p_v, filename = "TCGA-Omics_integration.png", dpi = 200, width = 20, height = 10)
} else {
  tcga_omics_sample <- read_delim(file = "TCGA_OMICS_INTER_SAMPLE_BARCODE.txt", delim = "\t") %>% 
    pull(1)
}


## exp
tcga_exp <- tcga_exp_raw %>%  
  dplyr::select(sample, all_of(tcga_omics_sample)) %>%
  mutate_if(is.numeric, function(value){
    normal <- 2^(value) - 0.001
    normal <- ifelse(normal < 0, 0, normal)
    log2(normal + 1) %>% return()
  }) %>% 
  separate(col = sample, into = c("Gene", "version"), sep = "\\.") %>% 
  dplyr::select(-version)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl_id <- tcga_exp %>% 
  dplyr::select(Gene) %>% 
  getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
        values = .$Gene, mart= mart)

tcga_exp_convert <- left_join(x = ensembl_id, y = tcga_exp , by = c("ensembl_gene_id" = "Gene")) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  dplyr::rename(Gene = hgnc_symbol) %>%
  filter(Gene != "") %>% 
  as_tibble()

exp_sd <- apply(X = tcga_exp_convert[, -1], FUN = sd, MARGIN = 1, na.rm = TRUE) %>% tibble(exp_sd = .)
exp_mean <- apply(X = tcga_exp_convert[, -1], FUN = mean, MARGIN = 1, na.rm = TRUE) %>% tibble(exp_mean = .)
exp_sd_mean <- bind_cols(tcga_exp_convert %>% select_at(1), exp_mean, exp_sd)

tcga_exp_index <- exp_sd_mean %>% 
  filter(exp_mean > 1 & exp_sd > 1) %>% 
  distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(Gene, Mean = exp_mean) %>%
  arrange(Gene)

save(tcga_exp_convert, file = "RData/TCGA-PANCAN-EXPRESSION.RData")
save(tcga_exp_index, file = "RData/TCGA-PANCAN-EXPRESSION_index.RData")

## mut
tcga_mut <- tcga_mut_raw %>% 
  filter(Tumor_Sample_Barcode %in% tcga_omics_sample) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  summarise(cnt = n()) %>% 
  mutate(cnt = 1)

tcga_samples <- tcga_mut %>% 
  pull(Tumor_Sample_Barcode) %>% 
  unique()

mut_longtowide <- lapply(X = tcga_samples, FUN = function(id){
  longtowide <- tcga_mut %>% 
    filter(Tumor_Sample_Barcode == id) %>%
    pivot_wider(names_from = "Tumor_Sample_Barcode", values_from = "cnt")
}) %>% 
  purrr::reduce(.x = ., full_join, by = "Hugo_Symbol") %>% 
  replace(is.na(.), 0)

# mutation frequency
tcga_mut_convert <- mut_longtowide %>% 
  mutate(MUT_TRUE = rowSums( . == 1), MUT_FALSE = rowSums(. == 0)) %>% 
  mutate(MUT_FREQ = MUT_TRUE / (MUT_TRUE + MUT_FALSE)) %>% 
  dplyr::select(Hugo_Symbol, MUT_FREQ) 
tcga_mut_convert$MUT_MEDIAN <- apply(mut_longtowide[, -1],1,function(v) median(as.numeric(v), na.rm = T))

tcga_mut_index <- tcga_mut_convert %>%
  filter(MUT_FREQ >= 0.01) %>% 
  select(Gene = Hugo_Symbol, Median = MUT_MEDIAN)

save(mut_longtowide, file = "RData/TCGA-PANCAN_MUTATION.RData")
# load("RData/TCGA-PANCAN_MUTATION.RData")
save(tcga_mut_index, file = "RData/TCGA-PANCAN_MUTATION_index.RData")


# Methylation
tcga_meth_convert <- tcga_meth_raw %>% 
  dplyr::select(Probe, all_of(tcga_omics_sample)) %>% 
  na.omit()

tcga_meth_index <- tcga_meth_convert %>% 
  mutate_if(is.numeric, .funs = function(col){
    ifelse(col <= 0.3, 1, 0) %>% return()
  }) %>% 
  mutate(METH_TRUE = rowSums( . == 1, na.rm = TRUE), METH_FALSE = rowSums(. == 0, na.rm = TRUE)) %>% 
  mutate(METH_FREQ = METH_TRUE / (METH_TRUE + METH_FALSE)) %>% 
  select(Probe, METH_FREQ) %>% 
  bind_cols(., apply(X = tcga_meth_convert[, -1], FUN = mean, MARGIN = 1, na.rm = TRUE) %>% 
              tibble(Mean = .)) %>%
  filter(METH_FREQ < 0.9) %>% 
  select(Probe, Mean)

save(tcga_meth_convert, file = "RData/TCGA-PANCAN_METHYLATION.RData")
save(tcga_meth_index, file = "RData/TCGA-PANCAN_METHYLATION_index.RData")


# CNA
tcga_cna <- tcga_cna_raw %>% 
  filter(Tumor_Sample_Barcode %in% tcga_omics_sample)

#  prep CNA function
tcga_sample <- unique(tcga_cna$Tumor_Sample_Barcode)
if (sum(tolower(colnames(tcga_cna)) %in% c("chr", "chromosome")) == 1) {
  col_idx <- which(tolower(colnames(tcga_cna)) %in% c("chr", "chromosome"))
  tcga_cna[which(tolower(tcga_cna[, col_idx]) == "x"), col_idx] <- 23
  tcga_cna[which(tolower(tcga_cna[, col_idx]) == "y"), col_idx] <- 24
  colnames(tcga_cna)[col_idx] <- "Chromosome"
}

dir.create("RData/TCGA-PANCAN-CNV", showWarnings = FALSE)
cnv_list <- list.files('RData/TCGA-PANCAN-CNV', full.names = T)
if(length(cnv_list) >= 22){
  bigTable_list <- list()
  for(index in 1:length(cnv_list)){
    print(cnv_list[index])
    load(cnv_list[index])
    bigTable_list[[index]] <- sample_bigTable %>% 
      lapply(X = ., FUN = function(col){ 
        col %>% select(-CNA)}) %>% 
      bind_cols(sample_bigTable[[1]] %>% select(CNA), .)
  }
} else {
  BP_SIZE <- 10 ^ 5
  chr <- unique(tcga_cna[, col_idx]) %>% pull()
  chrom_length <- c(249260000, 243200000, 198030000, 191160000, 
                    180920000, 171120000, 159140000, 146370000, 141220000, 
                    135540000, 135010000, 133860000, 115170000, 107350000, 
                    102540000, 90360000, 81200000, 78080000, 59130000, 63030000, 
                    48130000, 51310000, 155280000, 59380000)  ## chromomse bp length
  chrom_bin <- ceiling(chrom_length/BP_SIZE)
  bigTable <- data.frame(matrix(data = 0, ncol = length(tcga_sample) +  1, nrow = sum(chrom_bin)), 
                         stringsAsFactors = FALSE)
  colnames(bigTable) <- c("CNA", tcga_sample)
  
  k = 1
  for (i in 1:length(chr)) {
    bin_start <- seq(0, chrom_bin[i] - 1, 1)
    bin_end <- seq(1, chrom_bin[i], 1)
    bigTable$CNA[k:(k + chrom_bin[i] - 1)] <- paste(paste0("chr", i), paste0(bin_start, "to", bin_end), "10m", sep = "_")
    k = chrom_bin[i] + k
  }
  # not valid chr
  bigTable <- bigTable %>% 
    filter(CNA != "0")
  
  bigTable_Start_End <- bigTable %>% select(CNA) %>% 
    separate(col = "CNA", into = c("Chr", "Start_End", "bin"), sep = "_", remove = FALSE) %>% 
    select(-bin) %>% 
    separate(col = "Start_End", into = c("Start", "End"), sep = "to") %>% 
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  
  bigTable_list <- list()
  for (chr_ in chr[1:length(chr)]) { ## Chromosome 별로
    bigTable_chr <- bigTable %>% 
      filter(str_detect(string = CNA, pattern = paste0("^chr", chr_, "_")))
    cna.filterTable <- bigTable_Start_End %>% 
      filter(Chr == paste0("chr",chr_))
    idx.chr <- which(tcga_cna$Chromosome == chr_)
    Table.chr <- tcga_cna[idx.chr, ]
    cna.length.chr <- (Table.chr$End - Table.chr$Start) + 1  ## Start - end
    
    print(paste0("Chromosome : ", chr_))
    
    sample_bigTable <- pbmcapply::pbmclapply(X = 2:ncol(bigTable_chr), FUN = function(j){ #ncol(bigTable_chr)
      idx.cell <- which(Table.chr$Tumor_Sample_Barcode == colnames(bigTable_chr)[j])
      idx.big <- which(colnames(bigTable_chr) == colnames(bigTable_chr)[j])
      cellTable.chr <- Table.chr[idx.cell, ]
      cna.length <- cna.length.chr[idx.cell]
      
      tmp_bigTable <- bigTable_chr %>% select_at(.vars = c(1,idx.big))
      gc()
      
      for (k in cna.filterTable$End) {
        end_matrix <- data.frame(matrix(data = 1, nrow = nrow(cellTable.chr)) * BP_SIZE * k, stringsAsFactors = FALSE)
        end_matrix$cellTable <- cellTable.chr$End
        start_matrix <- data.frame(matrix(data = 1, nrow = nrow(cellTable.chr)) * BP_SIZE * (k - 1) + 1, stringsAsFactors = FALSE)
        start_matrix$cellTable <- cellTable.chr$Start
        overlap.length <- (BP_SIZE) + cna.length - (apply(end_matrix, 1, max) - apply(start_matrix, 1, min) + 1)
        overlap.length[overlap.length < 0] <- 0
        tmp_bigTable[k, 2] <- round(sum(overlap.length * cellTable.chr$Segment_Mean)/BP_SIZE, digits = 4)
      }
      
      return(tmp_bigTable)
    }, mc.cores = 24)
    
    # join
    save(sample_bigTable, file = paste0('RData/TCGA-PANCAN-CNV/chr', chr_, ".RData"))
    
    print("sample bigTable join...")
    bigTable_list[[chr_]] <- sample_bigTable %>% 
      lapply(X = ., FUN = function(col){ 
        col %>% select(-CNA)}) %>% 
      bind_cols(sample_bigTable[[1]] %>% select(CNA), .)
    print("sample bigTable join... : Done!")
    print(paste0("Chromosome : ", chr_, " Done!"))
  }
}

# all sample CNA
tcga_cna_convert <- bigTable_list %>% 
  bind_rows() %>% 
  separate(col = CNA, into = LETTERS[1:3], sep = "_", remove = FALSE) %>% 
  select(-B,-C) %>% 
  mutate(CHR = as.numeric(str_remove_all(A, "chr"))) %>% 
  select(-A) %>% 
  arrange(CHR) %>% 
  select(-CHR)

# filter
bigTable_com_zeros <- tcga_cna_convert %>% 
  mutate(CNA_NON_ZEROS = rowSums(. != 0, na.rm = TRUE), CNA_ZEROS = rowSums(. == 0, na.rm = TRUE)) %>% 
  mutate(CNA_ZEROS_FREQ = CNA_ZEROS / (CNA_NON_ZEROS + CNA_ZEROS)) %>% 
  select(CNA_ZEROS_FREQ)
cv <- apply(X = tcga_cna_convert[,-1], MARGIN = 1, function(x)sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100) %>% 
  tibble(CV = .)
absmean <- apply(X = tcga_cna_convert[,-1], MARGIN = 1, function(x) mean(abs(x), na.rm = TRUE)) %>% 
  tibble(ABS_MEAN = .)

tcga_cna_index <- tcga_cna_convert %>% 
  select(CNA) %>% 
  bind_cols(., bigTable_com_zeros, cv, absmean) %>% 
  filter(CNA_ZEROS_FREQ <= 0.05 & CV > 0.2 & ABS_MEAN > 0.15)

tcga_cna_index <-  tcga_cnv_index %>% 
  select(CNA) %>% 
  separate(col = CNA, into = c("Chr", "Start_END", "REMOVE1"), sep = "_", remove = FALSE) %>% 
  select(-REMOVE1) %>% 
  separate(col = Start_END, into = c("Start", "End"), sep = "to") %>% 
  mutate(Chr = as.numeric(str_remove_all(Chr, "chr")), Start = as.numeric(Start), End = as.numeric(End))

save(tcga_cna_convert, file = "RData/TCGA-PANCAN_CNA.RData")
save(tcga_cna_index, file = "RData/TCGA-PANCAN_CNA_index.RData")


# DeepDEP input source CCLE ----
## exp
tbl(con, "CCLE_expression") %>% 
  collect() %>% 
  rename(Gene = GENE) %>% 
  copy_to(dest = con, df = ., name = "DeepDEP_Expression", overwrite = T, temporary = F, indexes = list("Gene"))

## mut - nonsense, missense, frameshit IS, DEL, splice
mutation_type <- c('Nonsense_Mutation', 'Missense_Mutation', 'Frame_Shift_Del', 
                   'Frame_Shift_Ins','Splice_Site' )
mut <- tbl(con, "CCLE_mutation") %>% 
  collect() %>% 
  filter(Variant_Classification %in% mutation_type) %>% 
  group_by(DepMap_ID, Hugo_Symbol) %>% 
  summarise(cnt = n()) %>% 
  mutate(cnt = 1)

ccls <- mut %>% 
  pull(DepMap_ID) %>% 
  unique()

mut_longtowide <- lapply(X = ccls, FUN = function(id){
  longtowide <- mut %>% 
    filter(DepMap_ID == id) %>%
    pivot_wider(names_from = "DepMap_ID", values_from = "cnt")
}) %>% 
  reduce(.x = ., full_join, by = "Hugo_Symbol") %>% 
  replace(is.na(.), 0)

mut_longtowide %>% 
  rename(Gene = Hugo_Symbol) %>% 
  copy_to(dest = con, df = ., name = "DeepDEP_Mutation", overwrite = T, temporary = F, indexes = list("Gene"))

## CNA
cna_seg <- read_csv(file = "DeepDEP/CCLE_22Q2/CCLE_segment_cn.csv") %>% 
  select(-Status, -Source) %>% 
  select(CCLE_name = DepMap_ID, Chromosome, Start, End, Num_Probes, Segment_Mean)

cna_seg %>% 
  copy_to(dest = con, df = ., name = "DeepDEP_CNA", overwrite = T, temporary = F, indexes = list("CCLE_name"))

# methylation
#increase file download timeout
options(timeout = 600)

if(!file.exists("CCLs_methylation_GSE68379.Rds")){
  #get raw data - idats, processed beta matrix, etc.
  getGEOSuppFiles("GSE68379")
  
  #list files
  idatFiles <- list.files("DeepDEP_input/CCLs_methylation_GSE68379", pattern = "idat.gz$", full = TRUE)
  
  #decompress individual idat files
  sapply(idatFiles, gunzip, overwrite = TRUE)
  
  #read idats and create RGSet
  RGSet <- read.metharray.exp("DeepDEP_input/CCLs_methylation_GSE68379")
  saveRDS(RGSet, "CCLs_methylation_GSE68379.Rds")
} else {
  RGSet <- readRDS("CCLs_methylation_GSE68379.Rds")
}

# phenotype
geoMat <- getGEO("GSE68379")
pD.all <- pData(geoMat[[1]])
col_name <- pD.all$characteristics_ch1.3

# get methylation beta
grSet <- preprocessQuantile(RGSet)
grBeta <- getBeta(grSet) %>% 
  as.data.frame()
colnames(grBeta) <- col_name
grBeta_NULL <- which(colnames(grBeta) == "", arr.ind = TRUE)
grBeta <- grBeta[, -grBeta_NULL]

grBeta %>% 
  rownames_to_column(var = "Probe") %>% 
  copy_to(dest = con, df = ., name = "DeepDEP_Methylation", overwrite = T, temporary = F, indexes = list("Probe"))

# Preprocessing CCLE for train ----
sample_info <- tbl(con, "sample_info") %>% 
  collect() %>% 
  select(DepMap_ID, stripped_cell_line_name, COSMICID) %>%
  filter(!is.na(COSMICID)) %>% 
  arrange(DepMap_ID)

# gene expression
ccle_exp <- tbl(con, "DeepDEP_Expression") %>% collect()
ccle_exp_col_name <- ccle_exp %>% 
  colnames() %>% 
  .[-1] %>% 
  tibble(DepMap_ID = .)
ccle_exp_col_name_com <- left_join(x = ccle_exp_col_name, y = sample_info, by = "DepMap_ID") %>% 
  pull(stripped_cell_line_name)
colnames(ccle_exp) <- c("Gene", ccle_exp_col_name_com)

# COSMIC ID exist
ccle_exp_cosmic <- which(!is.na(colnames(ccle_exp)), arr.ind = TRUE)
ccle_exp <- ccle_exp %>% 
  select_at(ccle_exp_cosmic) %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  as.data.frame()



# gene mutation
ccle_mut <- tbl(con, "DeepDEP_Mutation") %>% collect()
ccle_mut_col_name <- ccle_mut %>% 
  colnames() %>% 
  .[-1] %>% 
  tibble(DepMap_ID = .)
ccle_mut_col_name_com <- left_join(x = ccle_mut_col_name, y = sample_info, by = "DepMap_ID") %>% 
  pull(stripped_cell_line_name)
colnames(ccle_mut) <- c("Gene", ccle_mut_col_name_com)

# COSMIC ID exist
ccle_mut_cosmic <- which(!is.na(colnames(ccle_mut)), arr.ind = TRUE)
ccle_mut <- ccle_mut %>% select_at(ccle_mut_cosmic) %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  as.data.frame()



# copy number
ccle_cna <- tbl(con, "DeepDEP_CNA") %>% collect() %>% 
  left_join(x = ., y = sample_info, by = c("CCLE_name" = "DepMap_ID")) %>% 
  select(-CCLE_name, -COSMICID) %>% 
  select(CCLE_name = stripped_cell_line_name, everything()) %>% 
  filter(!is.na(CCLE_name)) %>% 
  as.data.frame()


# methylation
ccle_meth <- tbl(con, "DeepDEP_Methylation") %>% collect()
ccle_meth_col_name <- ccle_meth %>% colnames() %>% 
  lapply(X = ., FUN = function(value){
    str_split(string = value, pattern = " ") %>% 
      unlist() %>% 
      .[2] %>% 
      return()
  }) %>% unlist() %>% 
  tibble(COSMICID = .) %>% 
  mutate(COSMICID = as.integer(COSMICID)) %>% 
  .[-1, ]
ccle_meth_col_name_com <- left_join(x = ccle_meth_col_name, y = sample_info, by = "COSMICID") %>% 
  pull(stripped_cell_line_name)
colnames(ccle_meth) <- c("Probe", ccle_meth_col_name_com)

# COSMIC ID exist
ccle_meth_cosmic <- which(!is.na(colnames(ccle_meth)), arr.ind = TRUE)
ccle_meth <- ccle_meth %>% select_at(ccle_meth_cosmic) %>% 
  distinct(Probe, .keep_all = TRUE) %>% 
  as.data.frame()



# gene dependency
ccle_gene_dependency <- tbl(con, "CRISPR_gene_effect") %>% collect()
ccle_gene_dependency_col_name <- ccle_gene_dependency %>% 
  colnames() %>% 
  .[-1] %>% 
  tibble(DepMap_ID = .)
ccle_gene_dependency_col_name_com <- left_join(x = ccle_gene_dependency_col_name, y = sample_info, by = "DepMap_ID") %>% 
  pull(stripped_cell_line_name)
colnames(ccle_gene_dependency) <- c("Gene", ccle_gene_dependency_col_name_com)

# COSMIC ID exist
ccle_gene_dependency_cosmic <- which(!is.na(colnames(ccle_gene_dependency)), arr.ind = TRUE)
ccle_gene_dependency <- ccle_gene_dependency %>% select_at(ccle_gene_dependency_cosmic) %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  as.data.frame()



# CCLE omics integration venn diagram
omics_ccl_list <- list(
  GeneDependency = ccle_gene_dependency %>% colnames() %>% .[-1],
  Exp = ccle_exp %>% colnames() %>% .[-1],
  Mutation = ccle_mut %>% colnames() %>% .[-1],
  Cna = ccle_cna %>% pull(CCLE_name) %>% unique(),
  Methylation = ccle_meth %>% colnames() %>% .[-1])
ccle_omics_intersection <- purrr::reduce(omics_ccl_list, intersect)

p_v <- ggVennDiagram::ggVennDiagram(x = omics_ccl_list) +
  ggtitle("CCLE-Omics integration") + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(p_v, filename = "CCLE-Omics_integration.png", dpi = 200, width = 20, height = 10)

ccle_exp_com <- ccle_exp %>% select(Gene, any_of(ccle_omics_intersection))
ccle_mut_com <- ccle_mut %>% select(Gene, any_of(ccle_omics_intersection))
ccle_cna_com <- ccle_cna %>% filter(CCLE_name %in% ccle_omics_intersection)
ccle_meth_com <- ccle_meth %>% select(Probe, any_of(ccle_omics_intersection))
ccle_gene_dependency_com <- ccle_gene_dependency %>% select(Gene, any_of(ccle_omics_intersection))

# Save train dataset
save(ccle_exp_com, file = "Train/CCLE-COSMIC-EXPRESSION.RData")
save(ccle_mut_com, file = "Train/CCLE-COSMIC-MUTATION.RData")
save(ccle_cna_com, file = "Train/CCLE-COSMIC-CNA.RData")
save(ccle_meth_com, file = "Train/CCLE-COSMIC-METHYLATION.RData")
save(ccle_gene_dependency_com, file = "Train/CCLE-COSMIC-GENEDEPENDENCY.RData")





# Final, CCLE-TCGA intersection ----
# TCGA index RData load
load("RData/TCGA-PANCAN-EXPRESSION_index.RData")
load("RData/TCGA-PANCAN-EXPRESSION.RData")

load("RData/TCGA-PANCAN_MUTATION_index.RData")
load("RData/TCGA-PANCAN_MUTATION.RData")

load("RData/TCGA-PANCAN_CNA_index.RData")
load("RData/TCGA-PANCAN_CNA.RData")

load("RData/TCGA-PANCAN_METHYLATION_index.RData")
load("RData/TCGA-PANCAN_METHYLATION.RData")

# CCLE-COSMIC RData load
load("TRAIN/CCLE-COSMIC-EXPRESSION.RData")
load("TRAIN/CCLE-COSMIC-MUTATION.RData")
load("TRAIN/CCLE-COSMIC-CNA.RData")
load("TRAIN/CCLE-COSMIC-METHYLATION.RData")

#
gene_expression_gene <- list(CCLE_EXP = ccle_exp_com$Gene, TCGA_EXP = tcga_exp_index$Gene)
ccle_tcga_gene <- gene_expression_gene %>% purrr::reduce(., intersect) %>% unique()

mutation_gene <- list(CCLE_MUT = ccle_mut_com$Gene, TCGA_MUT = tcga_mut_index$Gene)
ccle_tcga_mut <- mutation_gene %>% purrr::reduce(., intersect) %>% unique()

meth_probe <- list(CCLE_METH = ccle_meth_com$Probe, TCGA_METH = tcga_meth_index$Probe)
ccle_tcga_meth <- meth_probe %>% purrr::reduce(., intersect) %>% unique()

# cna_index <- .PrepCNA_custom(cna.original = ccle_cna_com)

venn_diagram(gene_expression_gene, type = "EXP")
venn_diagram(mutation_gene, type = "MUT")
venn_diagram(meth_probe, type = "METH")

exp.index <- tcga_exp_index %>% 
  filter(Gene %in% ccle_tcga_gene)
mut.index <- tcga_mut_index %>% 
  filter(Gene %in% ccle_tcga_mut)
meth.index <- tcga_meth_index %>% 
  filter(Probe %in% ccle_tcga_meth)

# index
save(exp.index, file = "TCGA_INDEX/CUSTOM/ccle_exp_custom_5454.RData")
save(mut.index, file = "TCGA_INDEX/CUSTOM/ccle_mut_custom_4946.RData")
save(meth.index, file = "TCGA_INDEX/CUSTOM/ccle_meth_custom_6231.RData")

# txt file
tcga_exp_ccl_paired <- exp.index %>% 
  select(Gene) %>% 
  left_join(x = ., y = tcga_exp_convert, by = "Gene")
tcga_exp_ccl_paired %>% 
  data.table::fwrite(file = "TCGA_INDEX/CUSTOM/tcga_exp_data_paired_with_ccl_custom.txt", sep = "\t")

tcga_mut_ccl_paired <- mut.index %>% 
  select(Gene) %>% 
  left_join(x = ., y = mut_longtowide, by = c("Gene" = "Hugo_Symbol"))
tcga_mut_ccl_paired %>% 
  data.table::fwrite(file = "TCGA_INDEX/CUSTOM/tcga_mut_data_paired_with_ccl_custom.txt", sep = "\t")

tcga_cna_ccl_paired <- tcga_cna_index %>% 
  select(CNA) %>% 
  left_join(x = ., y = tcga_cna_convert, by = "CNA")
tcga_cna_ccl_paired %>% 
  data.table::fwrite(file = "TCGA_INDEX/CUSTOM/tcga_cna_data_paired_with_ccl_custom.txt", sep = "\t")

tcga_meth_paired <- meth.index %>% 
  select(Probe) %>% 
  left_join(x = ., y = tcga_meth_convert, by = "Probe")
tcga_meth_paired %>% 
  data.table::fwrite(file = "TCGA_INDEX/CUSTOM/tcga_meth_data_paired_with_ccl_custom.txt", sep = "\t")

# TCGA-CCLE venn diagram
gene_expression_gene <- list(CCLE_EXP = ccle_exp_com$Gene, TCGA_EXP = tcga_exp_index$Gene)
mutation_gene <- list(CCLE_MUT = ccle_mut_com$Gene, TCGA_MUT = tcga_mut_index$Gene)
meth_probe <- list(CCLE_METH = ccle_meth_com$Probe, TCGA_METH = tcga_meth_index$Probe)

venn_diagram(gene_expression_gene, type = "EXP")
venn_diagram(mutation_gene, type = "MUT")
venn_diagram(meth_probe, type = "METH")

