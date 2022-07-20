source("src/function.R")

# train data
load("TRAIN/CCLE-COSMIC-EXPRESSION.RData")
load("TRAIN/CCLE-COSMIC-MUTATION.RData")
load("TRAIN/CCLE-COSMIC-CNA.RData")
load("TRAIN/CCLE-COSMIC-METHYLATION.RData")
load("TRAIN/CCLE-COSMIC-GENEDEPENDENCY.RData")

# for prediction
test_col <- sample(colnames(ccle_exp_com), 10)
# test_col <- c("Gene", test_col) %>% unique()

# train
## custom
Prep4DeepDEP_custom(
  exp.data = ccle_exp_com,
  mut.data = ccle_mut_com,
  meth.data = ccle_meth_com,
  cna.data = ccle_cna_com,
  dep.data = ccle_gene_dependency_com,
  mode = "training",
  filename.out = "train"
  )

# predict
## custom
Prep4DeepDEP_custom(exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
                    mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
                    meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
                    cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
                    # dep.data = ccle_gene_dependency_com %>% select(Gene, all_of(test_col)),
                    mode = "prediction",
                    filename.out = "predict")

Prep4DeepDEP(
  exp.data = ccle_exp_com %>% select(Gene, all_of(test_col)),
  mut.data = ccle_mut_com %>% select(Gene, all_of(test_col)),
  meth.data = ccle_meth_com %>% select(Probe, all_of(test_col)),
  cna.data = ccle_cna_com %>% filter(CCLE_name %in% test_col),
  mode = "prediction",
  filename.out = "predict")


