source("src/function.R")

path <- "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/PREDICTION"
now_date <- Sys.time() %>% str_split(" ") %>% unlist() %>% .[1]
save_path <- paste0(path, "/", now_date)

dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

ccls_wmbio <- read_delim("https://s3.us-west-2.amazonaws.com/secure.notion-static.com/0ccd993c-0f2c-4529-925f-40669970c4ed/WMBIO_UNIQUE_CELL_LINE.txt?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20220725%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220725T064945Z&X-Amz-Expires=86400&X-Amz-Signature=cebc944049e5556bf73499e22309d9065fbbd73e61dddda244eb6e8447413be1&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22WMBIO_UNIQUE_CELL_LINE.txt%22&x-id=GetObject",
                         delim = "\t", col_names = F) %>% 
  pull(1)


ccle_omics_extraction(ccls = cclw_wmbio,
                      CCLE_SAMPLE_INFO = "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/RAW/CCLs/sample_info.csv", 
                      CCLE_EXP_PATH = "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/RAW/CCLs/CCLE_expression.csv",
                      CCLE_MUT_PATH = "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/RAW/CCLs/CCLE_mutations.csv",
                      CCLE_METH_PATH = "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/RAW/CCLs/CCLs_methylation_GSE68379.Rds",
                      CCLE_CNA_PATH = "/home/wmbio/WORK/gitworking/DeepDEP/preprocessing/RAW/CCLs/CCLE_segment_cn.csv",
                      save_path = save_path)
