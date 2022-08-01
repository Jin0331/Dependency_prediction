source("src/function.R")

path <- "/home/wmbio/WORK/gitworking/Dependency_prediction//preprocessing/PREDICTION"
now_date <- Sys.time() %>% str_split(" ") %>% unlist() %>% .[1]
save_path <- paste0(path, "/", now_date)

dir.create(save_path, showWarnings = TRUE, recursive = TRUE)

ccls_wmbio <- read_delim("WMBIO_CCLS.txt", delim = "\t", col_names = F) %>% pull(1)

ccle_omics_extraction(ccls = cclw_wmbio,
                      CCLE_SAMPLE_INFO = "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/RAW/CCLs/sample_info.csv", 
                      CCLE_EXP_PATH = "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/RAW/CCLs/CCLE_expression.csv",
                      CCLE_MUT_PATH = "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/RAW/CCLs/CCLE_mutations.csv",
                      CCLE_METH_PATH = "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/RAW/CCLs/CCLs_methylation_GSE68379.Rds",
                      CCLE_CNA_PATH = "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/RAW/CCLs/CCLE_segment_cn.csv",
                      save_path = save_path)

# predict
Prep4DeepDEP(
  exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
  mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
  meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
  cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
  mode = "prediction",
  filename.out = paste0(save_path, "/predict_wmbio_ccls"))

