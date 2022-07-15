# CCLE omics integration venn diagram
venn_diagram <- function(inter_list, type = ""){
  p_v <- ggVennDiagram::ggVennDiagram(x = inter_list) +
    ggtitle("CCLE-TCGA omics integration") + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  ggsave(p_v, filename = paste0("CCLE-TCGA_integration_", type, ".png"), dpi = 200, width = 20, height = 10)
  return(p_v)
}

Prep4DeepDEP_custom <- function (exp.data = NULL, mut.data = NULL, meth.data = NULL, cna.data = NULL, dep.data = NULL, 
                          mode = c("Training", "Prediction"), filename.out = "data",
                          TCGA_INDEX_PATH = "TCGA_INDEX/"){
  
    check.cellNames <- NULL
    cat("Mode", mode, "\n\n")
    data_save_path <- Sys.time() %>% str_split(pattern = " ") %>% unlist() %>% .[1] %>% 
      paste0("DATA/", ., "/")
    dir.create(data_save_path, recursive = TRUE,showWarnings = FALSE)  
    
    # TCGA_INDEX_PATH <- "TCGA_INDEX/"
    
    if (sum(is.null(exp.data), is.null(mut.data), is.null(meth.data), is.null(cna.data)) == 4) {
      stop(c("All genomic profiles are missing. Please provide at least one of mut.data, exp.data, meth.data, and cna.data."))
    }
    if (is.null(dep.data) & tolower(mode) == "prediction") {
      cat("dep.data is not provided, running with the default 1298 DepOIs...", "\n")
      load(paste0(TCGA_INDEX_PATH, "CUSTOM/custom_dep_genes_1227.RData"))
    }
    if (is.null(dep.data) & tolower(mode) == "training") {
        cat("dep.data is not provided. Please provide gene dependency scores for the training mode...", "\n")
    }
    if (ncol(dep.data) == 1 & tolower(mode) == "training") {
        stop(c("Only one column detected in dep.data. Please provide gene dependency symbols and scores for the training mode."), call. = FALSE)
    }

    load(paste0(TCGA_INDEX_PATH, "CUSTOM/gene_fingerprints_CGP.RData")) # gene signature
    list.genes <- Prep4DeepDEP::.CheckGeneSymbol(dep.data = dep.data, filename.out = filename.out)
    n <- nrow(list.genes)
    
    # Gene expression
    if (!is.null(exp.data)) {
        if (!is.character(exp.data[1, 1])) {
            stop("exp.data format error, please check!", call. = FALSE)
        }
        cat(c("Exp started..."), "\n")
        colnames(exp.data)[1] <- "Gene"
        ncell <- ncol(exp.data[, -1])
        if (is.null(check.cellNames)) {
            check.cellNames <- colnames(exp.data[, -1])
        }else if (length(check.cellNames) != length(colnames(exp.data[, -1])) | sum(check.cellNames %in% colnames(exp.data[, -1])) != length(check.cellNames)) {
            stop(c("Cell line names are inconsistent!"), call. = FALSE)
        }
        cat("Precessing", paste(length(check.cellNames), "cell lines..."), "\n")
        inputData <- exp.data[!duplicated(exp.data$Gene), ]

        ######## TCGA - CCLE intersection        
        load(paste0(TCGA_INDEX_PATH, "CUSTOM/ccle_exp_custom_5454.RData"))

        # outputData <- merge(exp.index, inputData, by = "Gene", sort = FALSE, all.x = TRUE)
        outputData <- merge(exp.index, inputData, by = "Gene", sort = FALSE, all.x = TRUE)

        Gene <- outputData$Gene
        rownames(outputData) <- outputData$Gene
        value_NA <- rowSums(outputData[, -c(1, 2)])
        cat(sum(is.na(value_NA)), "genes with NA values in exp.data. Substitute by mean values of CCLE.", "\n")
        if (round((sum(is.na(value_NA))/nrow(outputData)), digits = 2) > 0.2) {
            warning("NA found in >20% genes, please check if input format is correct!")
        }
        for (i in 1:nrow(outputData)) {
            if (is.na(value_NA[i])) {
                outputData[i, is.na(outputData[i, ])] <- outputData$Mean[i]
            }
        }
        outputData <- round(as.matrix(outputData[, -1]), digits = 4)
        outputData.final.exp <- cbind(Gene, as.data.frame(outputData[, -1], stringsAsFactors = FALSE))
        outputData.final.exp <- outputData.final.exp[, c("Gene", check.cellNames)]
        if (tolower(mode) == "prediction") {
            data.table::fwrite(outputData.final.exp, 
                               file = paste0(data_save_path, paste(filename.out, "exp_prediction.txt", sep = "_")), sep = "\t", 
                               row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        if (tolower(mode) == "training") {
          k = 2:ncol(outputData.final.exp)
          rep_col_list <- pbmcapply::pbmclapply(X = k, FUN = function(index){
            rep_col.1 <- do.call("cbind", replicate(n, outputData.final.exp[, index], simplify = FALSE)) # 
            colnames(rep_col.1) <- paste0("C", index-1, "G", seq(1, n, 1))
            return(rep_col.1)
          }, mc.cores = 4)
          
          rep_col <- rep_col_list %>% 
            bind_cols(tibble(Gene = rownames(outputData.final.exp)), .) 
            
          data.table::fwrite(rep_col, file = paste0(data_save_path, paste(filename.out, "exp_training.txt", sep = "_")), 
                      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        cat("Exp completed!", "\n\n")
    }
    # Mutation
    if (!is.null(mut.data)) {
      if (!is.character(mut.data[1, 1])) {
        stop("mut.data format error, please check!", call. = FALSE)
      }
      cat(c("Mut started..."), "\n")
      colnames(mut.data)[1] <- "Gene"
      ncell <- ncol(mut.data[, -1])
      if (is.null(check.cellNames)) {
        check.cellNames <- colnames(mut.data[, -1])
      } else if (length(check.cellNames) != length(colnames(mut.data[, -1])) | sum(check.cellNames %in% colnames(mut.data[, -1])) != length(check.cellNames)) {
        stop(c("Cell line names are inconsistent!"), call. = FALSE)
      }
      cat("Precessing", paste(length(check.cellNames), "cell lines..."), "\n")
      inputData <- mut.data[!duplicated(mut.data$Gene), ]
      
      ######## TCGA - CCLE intersection        
      load(paste0(TCGA_INDEX_PATH, "CUSTOM/ccle_mut_custom_4946.RData"))
      
      outputData <- merge(mut.index, inputData, by = "Gene", sort = FALSE, all.x = TRUE)
      Gene <- outputData$Gene
      rownames(outputData) <- outputData$Gene
      value_NA <- rowSums(outputData[, -c(1, 2)])
      cat(sum(is.na(value_NA)), "genes with NA values in mut.data. Substitute by median values of CCLE.", "\n")
      if (round((sum(is.na(value_NA))/nrow(outputData)), digits = 2) > 0.2) {
        warning("NA found in >20% genes, please check if input format is correct!")
      }
      for (i in 1:nrow(outputData)) {
        if (is.na(value_NA[i])) {
          outputData[i, is.na(outputData[i, ])] <- outputData$Median[i]
        }
      }
      outputData.final.mut <- cbind(Gene, as.data.frame(outputData[, -c(1, 2)]))
      outputData.final.mut <- outputData.final.mut[, c("Gene", check.cellNames)]
      if (tolower(mode) == "prediction") {
        data.table::fwrite(outputData.final.mut, 
                           file = paste0(data_save_path, paste(filename.out, "mut_prediction.txt", sep = "_")), 
                           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
      if (tolower(mode) == "training") {
        k = 2:ncol(outputData.final.mut)
        rep_col_list <- pbmcapply::pbmclapply(X = k, FUN = function(index){
          rep_col.1 <- do.call("cbind", replicate(n, outputData.final.mut[, index], simplify = FALSE)) # 
          colnames(rep_col.1) <- paste0("C", index-1, "G", seq(1, n, 1))
          return(rep_col.1)
        }, mc.cores = 4)
        
        rep_col <- rep_col_list %>% 
          bind_cols(tibble(Gene = rownames(outputData.final.mut)), .) 
        
        data.table::fwrite(rep_col, 
                           file = paste0(data_save_path,paste(filename.out, "mut_training.txt", sep = "_")), 
                           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
      cat("Mut completed!", "\n\n")
    }
    # Methylation
    if (!is.null(meth.data)) {
        if (!is.character(meth.data[1, 1])) {
            stop("meth.data format error, please check!", call. = FALSE)
        }
        cat(c("Meth started..."), "\n")
        colnames(meth.data)[1] <- "Probe"
        ncell <- ncol(meth.data[, -1])
        if (is.null(check.cellNames)) {
            check.cellNames <- colnames(meth.data[, -1])
        } else if (length(check.cellNames) != length(colnames(meth.data[, -1])) | sum(check.cellNames %in% colnames(meth.data[, -1])) != length(check.cellNames)) {
            stop(c("Cell line names are inconsistent!"), call. = FALSE)
        }
        
        cat("Precessing", paste(length(check.cellNames), "cell lines..."), "\n")
        inputData <- meth.data[!duplicated(meth.data$Probe), ]
        
        load(paste0(TCGA_INDEX_PATH, "CUSTOM/ccle_meth_custom_6231.RData"))
        
        outputData <- merge(meth.index, inputData, by = "Probe", sort = FALSE, all.x = TRUE)
        Probe <- outputData$Probe
        rownames(outputData) <- outputData$Probe
        value_NA <- rowSums(outputData[, -c(1, 2)])
        cat(sum(is.na(value_NA)), "genes with NA values in meth.data. Substitute by 0.",  "\n")
        if (round((sum(is.na(value_NA))/nrow(outputData)), digits = 2) > 0.2) {
            warning("NA found in >20% genes, please check if input format is correct!")
        }
        for (i in 1:nrow(outputData)) {
            if (is.na(value_NA[i])) {
              if (sum(is.na(outputData[i, ])) == sum(ncol(outputData) - 2)) {
                  outputData[i, is.na(outputData[i, ])] <- 0
                } else {
                  outputData[i, is.na(outputData[i, ])] <- 0
                }
            }
        }
        outputData <- round(as.matrix(outputData[, -1]), digits = 4)
        outputData.final.meth <- cbind(Probe, as.data.frame(outputData[, -1]))
        outputData.final.meth <- outputData.final.meth[, c("Probe", check.cellNames)]
        if (tolower(mode) == "prediction") {
            data.table::fwrite(outputData.final.meth, 
                               file = paste0(data_save_path, paste(filename.out, "meth_prediction.txt", sep = "_")), 
                               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        if (tolower(mode) == "training") {
          k = 2:ncol(outputData.final.meth)
          rep_col_list <- pbmcapply::pbmclapply(X = k, FUN = function(index){
            rep_col.1 <- do.call("cbind", replicate(n, outputData.final.meth[, index], simplify = FALSE)) # 
            colnames(rep_col.1) <- paste0("C", index-1, "G", seq(1, n, 1))
            return(rep_col.1)
          }, mc.cores = 4)
          
          rep_col <- rep_col_list %>% 
            bind_cols(tibble(Probe = rownames(outputData.final.meth)), .) 
          
          data.table::fwrite(rep_col, 
                             file = paste0(data_save_path, paste(filename.out, "meth_training.txt", sep = "_")), 
                             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        cat("Meth completed!", "\n\n")
    }
    # Gene dependency
    if (!is.character(dep.data[1, 1])) {
        stop("dep.data format error, please check!", call. = FALSE)
    } else {
        cat(c("Fingerprint started..."), "\n")
        colnames(dep.data)[1] <- "Gene"
        idx <- which(fingerprint[1, ] %in% c("GeneSet", list.genes$Gene))
        outputData <- fingerprint[, idx]
        if (tolower(mode) == "prediction") {
            data.table::fwrite(outputData, 
                               file = paste0(data_save_path, paste(filename.out, "fingerprint_prediction.txt", sep = "_")), 
                               sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        if (tolower(mode) == "training") {
            outputData.train <- cbind(outputData[, 1], 
                                      do.call("cbind", replicate(ncell, outputData[, -1], simplify = FALSE)))
            outputData.train[1, ] <- c("GeneSet", paste0(paste0("C", rep(seq(1, ncell, 1), each = n)), "G", seq(1, n, 1)))
            data.table::fwrite(outputData.train, 
                               file = paste0(data_save_path, paste(filename.out, "fingerprint_training.txt", sep = "_")), 
                               sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
            
        }
        cat("Fingerprint completed!", "\n\n")
    }
    if (!is.null(dep.data) & tolower(mode) == "training") {
        if (is.null(check.cellNames)) {
            check.cellNames <- colnames(dep.data[, -1])
        }else if (length(check.cellNames) != length(colnames(dep.data[, -1])) | sum(check.cellNames %in% colnames(dep.data[, -1])) != length(check.cellNames)) {
            stop(c("Cell line names are inconsistent!"), call. = FALSE)
        }
        cat("Gene dependency scores (training mode) start...", "\n")
        crispr.input <- dep.data[which(dep.data$Gene %in% list.genes$Gene), 
                                 which(colnames(dep.data) %in% c("Gene", check.cellNames))]
        k = 2:ncol(crispr.input)
        crispr.output <- pbmcapply::pbmclapply(X = k, FUN = function(index){
          table <- as.data.frame(t(crispr.input[, index]))
          colnames(table) <- paste0("C", index - 1, "G", seq(1, n, 1))
          return(table)
        }, mc.cores = 3) %>% 
          bind_cols()
          
        crispr.output <- bind_cols(tibble(Dep_Score = "score"), crispr.output)
        data.table::fwrite(crispr.output, 
                           file = paste0(data_save_path, paste(filename.out, "DepScore_training.txt", sep = "_")), 
                           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        cat("Gene dependency scores (training) completed!", "\n\n")
    }
    if (!is.null(cna.data)) {
        if (sum(colnames(cna.data) %in% c("CCLE_name", "Chromosome", "Start", "End", "Num_Probe", "Segment_Mean")) != 5) {
            stop("cna.data format error, please check!", call. = FALSE)
        }
        cat(c("CNA started..."), "\n")
        ncell <- length(unique(cna.data$CCLE_name))
        if (is.null(check.cellNames)) {
            outputData.cna <- .PrepCNA_custom(cna.original = cna.data, filename.out, exportTable = FALSE)
        } else {
            idx <- which(cna.data$CCLE_name %in% check.cellNames)
            if (length(check.cellNames) != length(unique(cna.data$CCLE_name[idx])) | sum(check.cellNames %in% unique(cna.data$CCLE_name[idx])) != 
                  length(check.cellNames)) {
                stop(c("Cell line names are inconsistent!"), 
                  call. = FALSE)
            }
            outputData.cna <- .PrepCNA_custom(cna.original = cna.data[idx, ], filename.out, exportTable = FALSE)
            outputData.cna <- outputData.cna[, c("CNA", check.cellNames)]
        }
        if (tolower(mode) == "prediction") {
            colnames(outputData.cna)[1] <- "Bin"
            data.table::fwrite(outputData.cna, 
                               file = paste0(data_save_path, paste(filename.out, "cna_prediction.txt", sep = "_")), 
                               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        if (tolower(mode) == "training") {
          k = 2:ncol(outputData.cna)
          rep_col_list <- pbmcapply::pbmclapply(X = k, FUN = function(index){
            rep_col.1 <- do.call("cbind", replicate(n, outputData.cna[, index], simplify = FALSE)) # 
            colnames(rep_col.1) <- paste0("C", index-1, "G", seq(1, n, 1))
            return(rep_col.1)
          }, mc.cores = 4)
          
          rep_col <- rep_col_list %>% 
            bind_cols(tibble(Bin = outputData.cna$CNA), .) 
          
          data.table::fwrite(rep_col, 
                             file = paste0(data_save_path, paste(filename.out, "cna_training.txt", sep = "_")), 
                             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
        cat("CNA completed!", "\n\n")
    }
}

.PrepCNA_custom <- function (cna.original, filenames, filtered_table_path = "TCGA_INDEX/CUSTOM/cna_table_5915.RData", exportTable = FALSE) {
  library(progress)
  BP_SIZE <- 10 ^ 5
  cellLine <- unique(cna.original$CCLE_name)
  
  # chr type change
  cna.original <- cna.original %>%  
    mutate(Chromosome = ifelse(tolower(Chromosome) == "x", "23", 
                               ifelse(tolower(Chromosome) == "y", "24", Chromosome))) %>% 
    mutate(Chromosome = as.numeric(Chromosome))
  
  chr <- cna.original %>% pull(Chromosome) %>% unique()
  chrom_length <- c(249260000, 243200000, 198030000, 191160000, 
                    180920000, 171120000, 159140000, 146370000, 141220000, 
                    135540000, 135010000, 133860000, 115170000, 107350000, 
                    102540000, 90360000, 81200000, 78080000, 59130000, 63030000, 
                    48130000, 51310000, 155280000, 59380000)
  chrom_bin <- ceiling(chrom_length/BP_SIZE) # 10kb -> 100kb
  bigTable <- data.frame(matrix(data = 0, ncol = length(cellLine) + 1, nrow = sum(chrom_bin)), stringsAsFactors = FALSE)
  colnames(bigTable) <- c("CNA", cellLine)
  
  k = 1
  for (i in 1:length(chr)) {
    bin_start <- seq(0, chrom_bin[i] - 1, 1)
    bin_end <- seq(1, chrom_bin[i], 1)
    bigTable$CNA[k:(k + chrom_bin[i] - 1)] <- paste(paste0("chr", i), paste0(bin_start, "to", bin_end), "1m", sep = "_")
    k = chrom_bin[i] + k
  }
  t1 <- proc.time()
  load(filtered_table_path)
  
  for (i in unique(tcga_cnv_index$Chr)) {
    pb <- progress_bar$new(format = " Progress: [:bar] :percent, Estimated completion time: :eta",
                           total = ncol(bigTable), # totalnumber of ticks to complete (default 100)
                           clear = FALSE, # whether to clear the progress bar on completion (default TRUE)
                           width= 80 # width of the progress bar
    )
    print(paste0("Chromosome : ", i))
    cna.filterTable <- tcga_cnv_index[which(tcga_cnv_index$Chr == i), ]
    idx.chr <- which(cna.original$Chromosome == i)
    Table.chr <- cna.original[idx.chr, ]
    cna.length.chr <- (Table.chr$End - Table.chr$Start) + 1
    if (i == 1) {
      l = 0
    }
    else {
      l <- sum(chrom_bin[1:(i - 1)])
    }
    for (j in 2:ncol(bigTable)) {
      pb$tick()
      idx.cell <- which(Table.chr$CCLE_name == colnames(bigTable)[j])
      cellTable.chr <- Table.chr[idx.cell, ]
      cna.length <- cna.length.chr[idx.cell]
      for (k in cna.filterTable$End) {
        end_matrix <- data.frame(matrix(data = 1, nrow = nrow(cellTable.chr)) * BP_SIZE * k, stringsAsFactors = FALSE)
        end_matrix$cellTable <- cellTable.chr$End
        start_matrix <- data.frame(matrix(data = 1, nrow = nrow(cellTable.chr)) * BP_SIZE * (k - 1) + 1, stringsAsFactors = FALSE)
        start_matrix$cellTable <- cellTable.chr$Start
        overlap.length <- (BP_SIZE) + cna.length - (apply(end_matrix, 1, max) - apply(start_matrix, 1, min) + 1)
        overlap.length[overlap.length < 0] <- 0
        bigTable[l + k, j] <- round(sum(overlap.length * cellTable.chr$Segment_Mean)/BP_SIZE, digits = 4)
      }
    }
  }
  bigTable.filter <- merge(tcga_cnv_index, bigTable, by = "CNA", all.x = TRUE, sort = FALSE)
  bigTable <- bigTable.filter[, -c(2:4)]
  if (exportTable == TRUE) {
    data.table::fwrite(bigTable, file = paste(filenames, nrow(bigTable.filter), "_CNA_filter.txt", sep = "_"), sep = "\t", col.names = TRUE, 
                row.names = FALSE, quote = FALSE)
  }
  t2 <- proc.time()
  t <- round((t2 - t1)/60, digits = 2)
  print(c("Computation time (mins)", t[1:3]))
  return(bigTable)
}

