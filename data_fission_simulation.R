library(scDesign3)
library(tidyverse)
library(SingleCellExperiment)
library(parallel)
library(countsplit)
library(Seurat)
library(MASS)  # for rnegbin function
library(mclust)
#### seuart_processing function ####
seurat_processing <- function(countmat,
                              prj,
                              ncluster = 2,
                              resolution = 0.5,
                              FirstNClus = TRUE,
                              SCT = FALSE,
                              variable_feature = TRUE) {
  srtobj <- CreateSeuratObject(counts = countmat, project = prj, min.cells = 1, min.features = 1)
  all.genes <- rownames(srtobj)
  
  if(SCT) {
    srtobj <- SCTransform(srtobj, vst.flavor = "v2", verbose = FALSE)
  } else {
    srtobj <- NormalizeData(srtobj, normalization.method = "LogNormalize", scale.factor = 10000)
    srtobj <- FindVariableFeatures(srtobj, selection.method = "vst", nfeatures = 1000)
    srtobj <- ScaleData(srtobj, features = all.genes)
  }
  
  if(variable_feature) {
    srtobj <- RunPCA(srtobj, features = VariableFeatures(object = srtobj))
  } else {
    srtobj <- RunPCA(srtobj, features = all.genes)
  }
  srtobj <- FindNeighbors(srtobj, dims = 1:30)
  if (FirstNClus) {
    resolu <- 0.5
    srtobj <- FindClusters(srtobj, resolution = resolu)
    
    flagtime = 0
    while (length(unique(Idents(srtobj))) != ncluster & flagtime < 20) {
      if (length(unique(Idents(srtobj))) < ncluster & flagtime < 10) {
        resolu = resolu + 0.1
      } else if (length(unique(Idents(srtobj))) < ncluster & flagtime >= 10) {
        resolu = resolu + 0.01
      } else if (length(unique(Idents(srtobj))) > ncluster & flagtime < 3) {
        resolu = resolu - 0.10
      } else if (length(unique(Idents(srtobj))) > ncluster & flagtime >= 3) {
        resolu = resolu - 0.01
      }
      print(resolu)
      srtobj <- FindClusters(srtobj, resolution = resolu)
      flagtime = flagtime + 1
    }
    if (length(unique(Idents(srtobj))) == ncluster) {
      message(paste0("Successfully get ", ncluster, " clusters!"))
    } else {
      stop(paste0("Cannot make ", ncluster, " clusters!"))
    }
  } else {
    srtobj <- FindClusters(srtobj, resolution = resolution)
  }
  
  return(srtobj)
}
set.seed(123)


#### CalPvalueSimu function ####
CalPvalueSimu <- function(mat_exp, 
                          sigma_vec, 
                          true_label,
                          test = "wilcox",
                          variable_feature = TRUE, 
                          fastVersion = TRUE,
                          ifSparse = FALSE,
                          index = 1) {
  
  p_tbl <- tibble(Gene = rownames(mat_exp))
  
  ## Naive
  seurat.obj <- seurat_processing(countmat = mat_exp, 
                                  prj = "seurat",
                                  ncluster = 2,
                                  resolution = 0.5, 
                                  FirstNClus = TRUE, 
                                  SCT = FALSE, 
                                  variable_feature = variable_feature)
  label_naive = Idents(seurat.obj)
  ari_naive = adjustedRandIndex(true_label, label_naive)
  de_Seurat <- FindMarkers(seurat.obj, 
                           ident.1 = 0, 
                           ident.2 = 1, 
                           test.use = test,
                           logfc.threshold = 0, 
                           min.pct = 0, 
                           min.cells.feature = 0, 
                           min.cells.group = 0
  )
  
  p_Naive <- de_Seurat$p_val
  names(p_Naive) <- rownames(de_Seurat)
  p_Naive_tbl <- as_tibble(p_Naive, rownames = "Gene")
  p_Naive_tbl <- p_Naive_tbl %>% dplyr::rename(p_Naive = "value")
  print("Naive done!")
  ### CountSplit--data thinning
  overdisps.est <- sctransform::vst(mat_exp, n_genes = NULL)$model_pars[,1]
  input_data <- t(mat_exp)
  split <- countsplit::countsplit(X = input_data,folds = 2, epsilon = c(0.5, 0.5), overdisps = 1/sigma_vec)
  
  ##"train" "test"
  Xtrain <- t(split[[1]])
  Xtest <- t(split[[2]])
  
  exp_srtOj <- seurat_processing(Xtrain, prj = "train", FirstNClus = TRUE, ncluster = 2, SCT = FALSE, variable_feature = variable_feature)
  n_srtOj <- seurat_processing(Xtest, prj = "test", FirstNClus = TRUE, ncluster = 2, SCT = FALSE, variable_feature = variable_feature)
  
  ## Reset clusters
  Idents(n_srtOj) <- Idents(exp_srtOj)
  de_CountSplit <- FindMarkers(n_srtOj, 
                               ident.1 = 0, 
                               ident.2 = 1, 
                               test.use = test,
                               logfc.threshold = 0, 
                               min.pct = 0, 
                               min.cells.feature = 0, 
                               min.cells.group = 0
  )
  
  p_CountSplit <- de_CountSplit$p_val
  names(p_CountSplit) <- rownames(de_CountSplit)
  p_CountSplit_tbl <- as_tibble(p_CountSplit, rownames = "Gene")
  p_CountSplit_tbl <- p_CountSplit_tbl %>% dplyr::rename(p_CountSplit = "value")
  label_CountSplit = Idents(exp_srtOj)
  ari_CountSplit = adjustedRandIndex(true_label, label_CountSplit)
  print("Countsplit done!")
  
  ##### data fission ####
  simulate_process_givenX <- function(n, X, p = 0.5) {
    Z <- rbinom(n, X, p)
    return(Z)
  }
  n = dim(mat_exp)[2]
  
  Xtrain <- t(apply(mat_exp, 1, simulate_process_givenX, n = n))
  Xtest = mat_exp - Xtrain
  
  exp_srtOj <- seurat_processing(Xtrain, prj = "train", FirstNClus = TRUE, ncluster = 2, SCT = FALSE, variable_feature = variable_feature)
  n_srtOj <- seurat_processing(Xtest, prj = "test", FirstNClus = TRUE, ncluster = 2, SCT = FALSE, variable_feature = variable_feature)
  
  ## Reset clusters
  Idents(n_srtOj) <- Idents(exp_srtOj)
  de_data_fission <- FindMarkers(n_srtOj, 
                                 ident.1 = 0, 
                                 ident.2 = 1, 
                                 test.use = test,
                                 logfc.threshold = 0, 
                                 min.pct = 0, 
                                 min.cells.feature = 0, 
                                 min.cells.group = 0
  )
  
  p_data_fission <- de_data_fission$p_val
  names(p_data_fission) <- rownames(de_data_fission)
  p_data_fission_tbl <- as_tibble(p_data_fission, rownames = "Gene")
  p_data_fission_tbl <- p_data_fission_tbl %>% dplyr::rename(p_data_fission = "value")
  label_data_fission = Idents(exp_srtOj)
  ari_data_fission = adjustedRandIndex(true_label, label_data_fission)
  print("Data fission done!")
  
  ari_vec <- c(ari_naive, ari_CountSplit, ari_data_fission)
  
  pvalue <- Reduce(dplyr::left_join, list(p_Naive_tbl, p_CountSplit_tbl, p_data_fission_tbl))
  return(list(pvalue, ari_vec))
}


#### Set parameter ####
mu_para <- c(1, 0.25)
sigma_para <- c(1.5, 6)
gene_num_vec <- c(1000, 2000, 3000, 4000, 5000, 10000)
cell_num <- 1000
gene_module_size <- 50

## Fix corr gene number
corr_gene_num <- 100

for (gene_num in gene_num_vec) {
  res <- mclapply(seq_len(40), function(x) {
    
    DEcell <- sample(seq_len(cell_num), size = round(cell_num/2))
    true_label <- rep(0, cell_num)
    true_label[DEcell] <- 1
    gene_group_num <- corr_gene_num %/% gene_module_size
    corr_gene_index <- sample(1:round(0.5*gene_num), size = corr_gene_num)
    
    DEgene <- sample(setdiff(1:round(0.5*gene_num), corr_gene_index), size = 50) # 50 DE genes
    
    gene_name <- paste0("Gene", seq_len(gene_num))
    cell_name <- paste0("Cell", seq_len(cell_num))
    
    mu_vec <- rgamma(n = gene_num, shape = mu_para[1], rate = mu_para[2])
    mu_vec <- sort(mu_vec, decreasing = TRUE)
    sigma_vec <- rgamma(n = gene_num, shape = sigma_para[1], rate = sigma_para[2])
    
    corr_matrix_module <- lapply(1:gene_group_num, function(x) {
      toeplitz(c(1, seq(0.5, 0.0, length.out = gene_module_size-1)))})
    corr_matrix_sub <- Reduce(magic::adiag, corr_matrix_module)
    
    corr_matrix <- diag(nrow = gene_num)
    corr_matrix[corr_gene_index, corr_gene_index] <- corr_matrix_sub
    
    mean_mat <- matrix(mu_vec, nrow=cell_num, ncol=length(mu_vec), byrow=TRUE)
    for(j in DEgene) {
      mean_mat[-DEcell, j] <- (2^runif(1, -2, 2))*mean_mat[-DEcell, j]
    }
    
    
    colnames(mean_mat) <- gene_name
    rownames(mean_mat) <- cell_name
    
    sigma_mat <- matrix(sigma_vec, nrow=cell_num, ncol=length(sigma_vec), byrow=TRUE)
    
    colnames(sigma_mat) <- gene_name
    rownames(sigma_mat) <- cell_name
    
    rownames(corr_matrix) <- gene_name
    colnames(corr_matrix) <- gene_name
    zero_mat <- matrix(0, nrow = cell_num, ncol = gene_num)
    colnames(zero_mat) <- gene_name
    rownames(zero_mat) <- cell_name
    
    temp_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = t(mean_mat)))
    
    input_data <- data.frame(corr_group = as.factor(rep("1", cell_num)))
    rownames(input_data) <- cell_name
    
    colData(temp_sce) <- DataFrame(input_data)
    
    new_counts <- scDesign3::simu_new(sce = temp_sce, 
                                      assay_use = "counts", 
                                      mean_mat = mean_mat, 
                                      sigma_mat = sigma_mat, 
                                      #copula_list = list(NULL), 
                                      copula_list = list(`1` = corr_matrix),
                                      zero_mat = zero_mat, 
                                      n_cores = 1, 
                                      family_use = "nb", 
                                      new_covariate = NULL, 
                                      important_feature = "all",
                                      input_data = input_data,
                                      filtered_gene = NULL)
    
    mat_exp <- Matrix::Matrix(new_counts, sparse = TRUE)
    res = CalPvalueSimu(mat_exp,
                        sigma_vec = sigma_vec,
                        true_label = true_label,
                        index = x,
                        variable_feature = FALSE,
                        fastVersion = TRUE,
                        ifSparse = FALSE)
    res <- tryCatch(
      {CalPvalueSimu(mat_exp, 
                     sigma_vec = sigma_vec, 
                     index = x, 
                     true_label = true_label,
                     variable_feature = FALSE, 
                     fastVersion = TRUE, 
                     ifSparse = FALSE)},
      error = function(err) {
        return(NULL)
      }
    )
    res$DE <- FALSE
    res$Corr <- FALSE
    res$Corr[corr_gene_index] <- TRUE
    res$DE[DEgene] <- TRUE
    return(res)
  }, mc.cores = 20, mc.set.seed = TRUE, mc.preschedule = FALSE)
  saveRDS(res, file = paste0("data_path/comment", gene_num, ".rds"))
}

