library(tidyverse)
library(ggplot2)
library(wesanderson)
library(parallel)
theme_set(theme_bw())  # Set ggplot2 theme to black and white
palette.colors(palette = "Okabe-Ito")  # Set color palette
Vec_mtd <- c("Naive", "Data fission")
color_mtd <- c("blue", "red")
names(color_mtd) <- Vec_mtd

# Function to calculate False Discovery Proportion (FDP)
calFDP <- function(DisCover, DEgene) {
  DisCover <- as.numeric(stringr::str_remove(pattern = "Gene", DisCover))
  if(length(DisCover) == 0) {
    return(0)
  } else {
    length(setdiff(DisCover, DEgene))/length(DisCover)
  }
}

# Function to calculate Power
calPower <- function(DisCover, DEgene) {
  DisCover <- as.numeric(stringr::str_remove(pattern = "Gene", DisCover))
  if(length(DisCover) == 0) {
    return(0)
  } else {
    length(intersect(DisCover, DEgene))/length(DEgene)
  }
}
data_path = "data_path"  # Path to the data

# Function to create simulation plots
simulation_plot <- function(data_path) {
  n_cut <- 10
  FDR_threshold <- seq_len(n_cut)/100 # Define FDR thresholds
  
  common_pattern <- "comment"
  file_list <- list.files(data_path)
  # Extract gene numbers from file names
  gene_num_vec <- sapply(file_list, function(x) {
    x <- (gsub(pattern = common_pattern, replacement = "", x))
    x <- as.numeric(gsub(pattern = ".rds", replacement = "", x))
  })
  gene_num_vec <- sort(gene_num_vec)
  
  # Compute FDP for each gene number
  full_tbl_FDP <- mclapply(gene_num_vec, function(cutoff) {
    input_data <- readRDS(paste0(data_path, common_pattern, cutoff , ".rds"))
    input_data = compact(input_data)
    summary_list <- lapply(seq_len(length(input_data)), function(x) {
      example_tbl_list <- input_data[[x]]
      example_tbl = example_tbl_list[[1]]
      
      qvalue_tbl <- example_tbl %>% dplyr::mutate(q_Naive = p.adjust(p_Naive, method = "BH"),
                                                  q_data_fission = p.adjust(p_data_fission, method = "BH"))
      qvalue_tbl[is.na(qvalue_tbl)] <- 1 # Replace NA values with 1
      index_DE = which(example_tbl_list$DE==TRUE)
      qvalue_tbl$DE = FALSE
      qvalue_tbl$DE[index_DE] = TRUE
      # setting rownames for qvalue_tbl
      rownames(qvalue_tbl) = qvalue_tbl$Gene
      DEgene <- which((qvalue_tbl$DE))
      data_tbl <- suppressMessages(sapply(FDR_threshold, function(i) {
        FDP_DD <- calFDP(qvalue_tbl$Gene[which(qvalue_tbl$q_Naive <= i)], DEgene)
        FDP_data_fisson <- calFDP(qvalue_tbl$Gene[which(qvalue_tbl$q_data_fission <= i)], DEgene)
        res <- c(FDP_DD, FDP_data_fisson)
        names(res) <- c("Naive", "Data fission")
        res
      })) %>% t() %>% `rownames<-`(FDR_threshold)
      
      data_tbl
    })
    
    # Compute mean FDP
    mean_tbl <- Reduce("+", summary_list)/ length(summary_list)
    mean_tbl <- mean_tbl %>% as_tibble(rownames = "Nominal FDR") %>% tidyr::pivot_longer(-c("Nominal FDR"), names_to = "Method", values_to = "FDR") %>% 
      dplyr::mutate(Method = factor(Method, levels = c("Naive", "Data fission")), `Nominal FDR` = as.numeric(`Nominal FDR`)) %>% dplyr::mutate(Gene_num = cutoff, Type = "mean")
    
    res_tbl <- mean_tbl
    res_tbl
  }, mc.cores = length(gene_num_vec))
  
  full_tbl_Power <- mclapply(gene_num_vec, function(cutoff) {
    input_data <- readRDS(paste0(data_path, common_pattern, cutoff , ".rds"))
    input_data = compact(input_data)
    summary_list <- lapply(seq_len(length(input_data)), function(x) {
      example_tbl_list <- input_data[[x]]
      example_tbl = example_tbl_list[[1]]
      
      qvalue_tbl <- example_tbl %>% dplyr::mutate(q_Naive = p.adjust(p_Naive, method = "BH"),
                                                  q_data_fission = p.adjust(p_data_fission, method = "BH"))
      qvalue_tbl[is.na(qvalue_tbl)] <- 1
      index_DE = which(example_tbl_list$DE==TRUE)
      qvalue_tbl$DE = FALSE
      qvalue_tbl$DE[index_DE] = TRUE
      # setting rownames for qvalue_tbl
      rownames(qvalue_tbl) = qvalue_tbl$Gene
      DEgene <- which((qvalue_tbl$DE))
      data_tbl <- suppressMessages(sapply(FDR_threshold, function(i) {
        Power_DD <- calPower(qvalue_tbl$Gene[which(qvalue_tbl$q_Naive <= i)], DEgene)
        Power_data_fisson <- calPower(qvalue_tbl$Gene[which(qvalue_tbl$q_data_fission <= i)], DEgene)
        res <- c(Power_DD, Power_data_fisson)
        names(res) <- c("Naive", "Data fission")
        res
      })) %>% t() %>% `rownames<-`(FDR_threshold)

      data_tbl
    })
    
    mean_tbl <- Reduce("+", summary_list)/ length(summary_list)
    mean_tbl <- mean_tbl %>% as_tibble(rownames = "Nominal FDR") %>% tidyr::pivot_longer(-c("Nominal FDR"), names_to = "Method", values_to = "Power") %>% 
      dplyr::mutate(Method = factor(Method, levels = c("Naive", "Data fission")), `Nominal FDR` = as.numeric(`Nominal FDR`)) %>% dplyr::mutate(Gene_num = cutoff, Type = "mean")
    
    res_tbl <- mean_tbl
    res_tbl
  }, mc.cores = length(gene_num_vec))
  
  p_corr_prop_FDR <- full_tbl_FDP %>% bind_rows() %>% dplyr::mutate(FDR = 100*FDR, `Nominal FDR` = 100*`Nominal FDR`)%>% tidyr::pivot_wider(names_from = Type, values_from = FDR) %>% ggplot(aes(x = `Nominal FDR`, y = mean, color = Method)) + facet_wrap(~Gene_num, nrow = 1) + geom_line(aes(linetype = Method)) + geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 2)  + scale_y_continuous(breaks = seq_len(10)*10)+ scale_x_continuous(breaks = seq_len(n_cut)) +
    scale_linetype_manual(values=c(3, 4, 5, 6, 1)) + theme(aspect.ratio = 1) + xlab("Nominal FDR (%)") + ylab("Approximate FDR (%)") + scale_color_manual(values = color_mtd) + ylim(0, 100)
  
  
  p_corr_prop_Power <- full_tbl_Power %>% bind_rows() %>% dplyr::mutate(Power = 100*Power, `Nominal FDR` = 100*`Nominal FDR`)%>% tidyr::pivot_wider(names_from = Type, values_from = Power) %>% ggplot(aes(x = `Nominal FDR`, y = mean, color = Method)) + facet_wrap(~Gene_num, nrow = 1) + geom_line(aes(linetype = Method)) + scale_y_continuous(breaks = seq_len(10)*10)+ scale_x_continuous(breaks = seq_len(n_cut)) +
    scale_linetype_manual(values=c(3, 4, 5, 6, 1)) + theme(aspect.ratio = 1) + xlab("Nominal FDR (%)") + ylab("Power (%)") + scale_color_manual(values = color_mtd)+ ylim(0, 100)
  final_p <- ggpubr::ggarrange(plotlist = list(p_corr_prop_FDR+ ggtitle("Unsupervised setting"), p_corr_prop_Power), 
                               ncol = 1, 
                               common.legend = TRUE,
                               legend = "bottom")
  
  return(final_p)
}

# Save the final plot to a PDF file
pdf(paste0(data_path, "simulation_plot.pdf"), width = 10, height = 10)
simulation_plot(data_path)
dev.off()
