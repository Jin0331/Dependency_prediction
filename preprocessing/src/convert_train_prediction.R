source("~/WORK/gitworking/Dependency_prediction/preprocessing/src/function.R")

# train data
load("TRAIN/2022-07-27/CCLE-COSMIC-EXPRESSION.RData")
load("TRAIN/2022-07-27/CCLE-COSMIC-MUTATION.RData")
load("TRAIN/2022-07-27/CCLE-COSMIC-CNA.RData")
load("TRAIN/2022-07-27/CCLE-COSMIC-METHYLATION.RData")
load("TRAIN/2022-07-27/CCLE-COSMIC-GENEDEPENDENCY.RData")

# for prediction
test_col <- sample(colnames(ccle_exp_com), 3)
# test_col <- c("Gene", test_col) %>% unique()

# train
## custom
# Prep4DeepDEP_custom(
#   exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
#   mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
#   meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
#   cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
#   dep.data = ccle_gene_dependency_com %>% select(Gene, all_of(test_col)),
#   mode = "training",
#   filename.out = "train"
#   )

path <- "/home/wmbio/WORK/gitworking/Dependency_prediction/preprocessing/DATA"
now_date <- Sys.time() %>% str_split(" ") %>% unlist() %>% .[1]
save_path <- paste0(path, "/", now_date)
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# Prep4DeepDEP(
#     exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
#     mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
#     meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
#     cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
#     dep.data = ccle_gene_dependency_com %>% select(Gene, all_of(test_col)),
#   mode = "training",
#   filename.out = paste0(save_path, "/training_default")
# )

Prep4DeepDEP_custom(
  exp.data = ccle_exp_com ,
  mut.data = ccle_mut_com ,
  meth.data = ccle_meth_com ,
  cna.data = ccle_cna_com,
  dep.data = ccle_gene_dependency_com ,
  mode = "training",
  filename.out = paste0(save_path, "/training_custom")
)


# predict
## custom
# Prep4DeepDEP_custom(
#   exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
#   mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
#   meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
#   cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
#   # dep.data = ccle_gene_dependency_com %>% select(Gene, all_of(test_col)),
#   mode = "prediction",
#   filename.out = "predict"
# )

Prep4DeepDEP(
  exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
  mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
  meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
  cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
  mode = "prediction",
  filename.out = "predict")


