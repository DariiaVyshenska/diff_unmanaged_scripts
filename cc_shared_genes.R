getwd()
setwd("D:/GitHub/R_wd/")

library(gtools)

#======================
# FUNCTIONS
#======================


# FUNCTION: transforming a list into a data.frame for exprot:
# input - list of vectors with names;
# output - data.frame with names of vectors in column 1, collapsed comma sep.
# values from each vector in column 2
transform_list_out <- function(input_list, annotations){

  out_list <- data.frame(KD = character(), 
                         Genes=character(), 
                         stringsAsFactors=FALSE) 
  for(i in 1:length(input_list)){
    symbols <- unique(annotations$geneID[annotations$EnsemblID %in% input_list[[i]]])
    col_genes <- paste(symbols, collapse = ",")
    out_list[i,] <- c(names(input_list[i]), col_genes)
  }
  return(out_list)
}


# FUNCTION: importing data
# function takes a path to a folder that contains tab delimited files
# files are imported as data.frames and all of them put into a list 
# each table in the list has a name of the file name from wich it was imported
import_df_asList <- function(path){
  list_of_tables <- list()
  for (file in list.files(path)){
    #file <- "TPX2.txt"
    cat("Importing file: ", file, "\n", sep = "")
    file_path <- paste(input_folder, file, sep = "/")
    table <- read.table(file = file_path, 
                        sep = '\t', header = TRUE, check.names = F, stringsAsFactors = F) 
    gene_name <- gsub(".txt", "", file)
    list_of_tables[[gene_name]] <- table
  }
  return(list_of_tables)
}


# FUNCTION: filtering list of tables for PUC complience
# input - list of tables that contain column "PUC-complience?" (0 or 1)
# and filters this list to contain either all genes (puc = NA), 
# puc-complient genes (puc = T), puc not-complient genes (puc = F)
# output - list of tables
check_puc <- function(all_tables_test, puc_threshold = NA){
  if (is.na(puc_threshold)){
    cat("No filtering for PUC is done!\n")
    return(all_tables)
  } else if (puc_threshold == T){
    puc_val <- 1
  } else if (puc_threshold == F){
    puc_val <- 0
  } else if (cat("'puc' argument must be T/F/NA!\n"))
  cat("Filters PUC\n")
  new_tables_list <- list()
  for (table in 1:length(all_tables)){
    table_name <- names(all_tables[table])
    print(table_name)
    condition <- all_tables[[table]][["PUC-complient?"]] == puc_val
    new_tables_list[[table_name]] <- all_tables[[table]][eval(condition),]
  }
  return(new_tables_list)
}


# FUNCTION: filtering tables, creating list of genes and list of their ratios
# input - list of imported tables
# output - matrix with shared gene numbers/normalized ratios
#
# all_tables_lis - list of tables from BRB
# puc = NA/T/F
# pval_threshold = numeric OR F
# fdr_threshold = numeric OR F
# shared_direct = "same", "diff" or "all"
# normalized = F OR T
# export = T/F
# annotation_file_path = file to path containing two columns - EnsemblID and 
#geneID (tab del)
get_matrix <- function (all_tables_list, puc = NULL, pval_threshold, 
                        fdr_threshold, shared_direct, 
                        normalized, export = F, annotation_file_path = NA) {
  # filtering for PUC if nesessary
  if(!is.null(puc)){
    all_tables_list <- check_puc(all_tables = all_tables_list, puc_threshold = puc)
    }
  
  genes_list <- list()
  ratio_list <- list()
  
  # Extracting sublists based on p-value and FDR thresholds
  for(table in 1:length(all_tables_list)){
    sel_vec_signif <- vector()
    gene_name <- NA
    cat("Working with: ", names(all_tables_list[table]), "\n")
    if(is.numeric(pval_threshold) & is.numeric(fdr_threshold)){
      cat("Filtering by p-value and FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold &+
        all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (is.numeric(pval_threshold) & fdr_threshold == F){
      cat("Filtering only by p-value.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold
    } else if (pval_threshold == F & is.numeric(fdr_threshold)){
      cat("Filtering only by FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (pval_threshold == F & fdr_threshold == F){
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < 1.1
    } else { cat ("FDR and p-value MUST be either 'F' or numeric!\n")}
    
    gene_name <- names(all_tables_list[table])
    genes_list[[gene_name]] <- all_tables_list[[table]][["UniqueID"]][sel_vec_signif]
    ratio_list[[gene_name]] <- all_tables_list[[table]][["KD/C-all"]][sel_vec_signif]
  }
  
  cat(lengths(genes_list), "\n")
  
  ####
  # Writing out all genes that pass the threshold in an .csv file
  if(export == T) {
    if(!is.na(annotation_file_path)){
      annotation_file <- annotation_file_path
      annotations <- read.table(annotation_file_path, stringsAsFactors = F, 
                                check.names = F, header = T)
    } else {cat("You must provide annotation file!\n")
      break
      }
    

    genes_list_tr <- transform_list_out(genes_list, annotations)
    fullList_file_name <- paste("fullGeneListFDR-", 
                                fdr_threshold, "pval-", pval_threshold,".csv", sep="")
    write.csv(genes_list_tr, fullList_file_name, row.names = F)
  }
  #####
  
  # vector with total numbers of genes controlled by each KD - for normalization
  all_numbers <- lengths(genes_list)
  # collecting KD genes names
  nm_vec <- names(genes_list)
  # creating matrix to populate with numbers of shared genes
  shared_numbers <- matrix(0, nrow = length(genes_list), ncol = length(genes_list))
  rownames(shared_numbers) <- nm_vec
  colnames(shared_numbers) <- nm_vec
  
  # creating a list to populate with the shared genes names - for export in .csv file
  shared_genes <- list()
  
  # populating matrix and list based on directionality and normalization parameters
  for (i in 1:(length(genes_list) - 1)){
#i <- 1
    cat("I'm working with: ", nm_vec[i], "\n", sep = "")
    for (j in (i+1):length(genes_list)){
#j <- 2
      cat("comparing with ", names(genes_list[j]), "\n")
      colhead <- paste(names(genes_list[i]),"vs", names(genes_list[j]), sep = "")
      shared_only <- genes_list[[i]][genes_list[[i]] %in% genes_list[[j]]]
      length(shared_only)
      
      if(shared_direct == "all"){
        if (normalized == F){shared_numbers[i,j] <- length(shared_only)
        } else if (normalized == T){
          shared_numbers[i,j] <- length(shared_only)/mean(c(all_numbers[i], all_numbers[j]))
        } else {cat("ERROR: parameter 'normalized' must be T or F!\n")}
        if(length(shared_only) == 0){
          shared_genes[[colhead]] <- NA
        } else {shared_genes[[colhead]] <- shared_only}
      } else {
        shared_same_dir <- c()
        
        for(gene in shared_only){
          ratio_i <- ratio_list[[i]][genes_list[[i]] == gene]
          ratio_j <- ratio_list[[j]][genes_list[[j]] == gene]
          print(gene)
          print(ratio_i)
          print(ratio_j)

          if(shared_direct == "same"){
            condition <- expression(ratio_i > 1 & ratio_j > 1 | ratio_i < 1 & ratio_j < 1)
          } else if (shared_direct == "diff"){
            condition <- expression(ratio_i < 1 & ratio_j > 1 | ratio_i > 1 & ratio_j < 1)
          }
          #print(condition)
          #print(eval(condition))
          if(eval(condition) == T){
            shared_same_dir <- c(shared_same_dir, gene)
          } else if (eval(condition) == F){print("Different direction!")
            } else {print("ERROR in conditions!")}
        }
        
        #print(shared_same_dir)
        if (normalized == F){shared_numbers[i,j] <- length(shared_same_dir)
        } else if (normalized == T){
          shared_numbers[i,j] <- length(shared_same_dir)/mean(c(all_numbers[i], all_numbers[j]))
        } else {cat("ERROR: parameter 'normalized' must be T or F!\n")}
        if(length(shared_same_dir) == 0){
          shared_genes[[colhead]] <- NA
        } else {shared_genes[[colhead]] <- shared_same_dir}
      }
    } 
  }
  
  if(export == T) {
    # exporting in .csv files matrix with numbers/ratios of shared genes
    # and list of shared genes for each KD
  
    shared_g_df <- transform_list_out(shared_genes, annotations)
    matrix_file_name <- paste("shared-numb-fdr", 
                              fdr_threshold, "KDdir", toupper(shared_direct),".csv", sep="")
    gene_list_file <- paste("shared-genes-fdr", 
                            fdr_threshold, "KDdir", toupper(shared_direct),".csv", sep = "")
    write.csv(shared_numbers, matrix_file_name)
    write.csv(shared_g_df, gene_list_file, row.names = F)
  }
  return(shared_numbers)
}




# this function only gives out a matrix that was normalized in one of three 
# different ways. PUC, p-value and FDR filtering is performed within a function
# Normalizations:
# A - shared_same_dir/all_shared
# B - shared_same_dir/union_of_all_regulated
# C - (shared_same_dir - shared_different_dir)/union_of_all_regulated
# D - shared_different_dir/union_of_all_regulated
# E - all_shared/union_of_all_regulated


normaliz_matrix <- function (all_tables_list, puc = NULL, pval_threshold, 
                        fdr_threshold, normalization) {
#all_tables_list <- all_tables
#puc = NULL
#pval_threshold <- 0.001
#fdr_threshold <- F
#normalization <- "C"
  # filtering for PUC if nesessary
  if(!is.null(puc)){
    all_tables_list <- check_puc(all_tables = all_tables_list, puc_threshold = puc)
    }
  
  genes_list <- list()
  ratio_list <- list()
  
  # Extracting sublists based on p-value and FDR thresholds
  for(table in 1:length(all_tables_list)){
    sel_vec_signif <- vector()
    gene_name <- NA
    cat("Working with: ", names(all_tables_list[table]), "\n")
    if(is.numeric(pval_threshold) & is.numeric(fdr_threshold)){
      cat("Filtering by p-value and FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold &+
        all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (is.numeric(pval_threshold) & fdr_threshold == F){
      cat("Filtering only by p-value.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold
    } else if (pval_threshold == F & is.numeric(fdr_threshold)){
      cat("Filtering only by FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (pval_threshold == F & fdr_threshold == F){
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < 1.1
    } else { cat ("FDR and p-value MUST be either 'F' or numeric!\n")}
    
    gene_name <- names(all_tables_list[table])
    genes_list[[gene_name]] <- all_tables_list[[table]][["UniqueID"]][sel_vec_signif]
    ratio_list[[gene_name]] <- all_tables_list[[table]][["KD/C-all"]][sel_vec_signif]
  }
  
  cat(lengths(genes_list), "\n")
  # vector with total numbers of genes controlled by each KD - for normalization
  all_numbers <- lengths(genes_list)
  # collecting KD genes names
  nm_vec <- names(genes_list)
  # creating matrix to populate with numbers of shared genes
  shared_numbers <- matrix(0, nrow = length(genes_list), ncol = length(genes_list))
  rownames(shared_numbers) <- nm_vec
  colnames(shared_numbers) <- nm_vec
  

  # populating matrix and list based on directionality and normalization parameters
  for (i in 1:(length(genes_list) - 1)){
    #i <- 2
    cat("\nI'm working with: ", nm_vec[i], "\n", sep = "")
    for (j in (i+1):length(genes_list)){
      #j <- 2
      cat("comparing with ", names(genes_list[j]), "\n")
      colhead <- paste(names(genes_list[i]),"vs", names(genes_list[j]), sep = "")
      shared_only <- genes_list[[i]][genes_list[[i]] %in% genes_list[[j]]]
      all_shared <- length(shared_only)
      cat("Shared only is: ", shared_only, "\n")
      
      
        shared_same_dir <- c()
        shared_diff_dir <- c()
        for(gene in shared_only){
          ratio_i <- ratio_list[[i]][genes_list[[i]] == gene]
          ratio_j <- ratio_list[[j]][genes_list[[j]] == gene]
          print(gene)
          print(ratio_i)
          print(ratio_j)
      
          if(ratio_i > 1 & ratio_j > 1 | ratio_i < 1 & ratio_j < 1){
            # this two genes have same direction
            shared_same_dir <- c(shared_same_dir, gene)
          } else if (ratio_i < 1 & ratio_j > 1 | ratio_i > 1 & ratio_j < 1){
            # this two genes have different direction
            shared_diff_dir <- c(shared_diff_dir, gene)
          }
        }
          shared_same <- length(shared_same_dir)
          shared_diff <- length(shared_diff_dir)
          
          if(all_shared == 0){
            shared_numbers[i,j] <- 0
          } else(
          if(normalization == "A"){
            shared_numbers[i,j] <- shared_same/all_shared
            cat("shared number is:", shared_numbers[i,j], "\n")
            #shared_numbers[i,j] <- length(shared_same_dir)/mean(c(all_numbers[i], all_numbers[j]))
          } else if(normalization == "B"){
            shared_numbers[i,j] <- shared_same/(all_numbers[i]+all_numbers[j] - all_shared)
          } else if (normalization == "C"){
            shared_numbers[i,j] <- (shared_same - shared_diff)/(all_numbers[i] + all_numbers[j] - all_shared)
            cat("shared number is:", shared_numbers[i,j], "\n")
          } else if(normalization == "D"){
            shared_numbers[i,j] <- shared_diff/(all_numbers[i]+all_numbers[j] - all_shared)
          } else if(normalization == "E"){
            shared_numbers[i,j] <- all_shared/(all_numbers[i]+all_numbers[j] - all_shared)
          }else(print("You must specify normalization!"))
          )
    } 
  }

   return(shared_numbers)
}

# Function: same as normaliz_matrix but returns numbers, not ratios
#direction: same/diff/all/union
shared_numbers <- function (all_tables_list, puc = NULL, pval_threshold, 
                             fdr_threshold, direction) {
  #all_tables_list <- all_tables
  #puc = NULL
  #pval_threshold <- 0.001
  #fdr_threshold <- F
  #normalization <- "C"
  # filtering for PUC if nesessary
  if(!is.null(puc)){
    all_tables_list <- check_puc(all_tables = all_tables_list, puc_threshold = puc)
  }
  
  genes_list <- list()
  ratio_list <- list()
  
  # Extracting sublists based on p-value and FDR thresholds
  for(table in 1:length(all_tables_list)){
    sel_vec_signif <- vector()
    gene_name <- NA
    cat("Working with: ", names(all_tables_list[table]), "\n")
    if(is.numeric(pval_threshold) & is.numeric(fdr_threshold)){
      cat("Filtering by p-value and FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold &+
        all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (is.numeric(pval_threshold) & fdr_threshold == F){
      cat("Filtering only by p-value.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold
    } else if (pval_threshold == F & is.numeric(fdr_threshold)){
      cat("Filtering only by FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (pval_threshold == F & fdr_threshold == F){
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < 1.1
    } else { cat ("FDR and p-value MUST be either 'F' or numeric!\n")}
    
    gene_name <- names(all_tables_list[table])
    genes_list[[gene_name]] <- all_tables_list[[table]][["UniqueID"]][sel_vec_signif]
    ratio_list[[gene_name]] <- all_tables_list[[table]][["KD/C-all"]][sel_vec_signif]
  }
  
  cat(lengths(genes_list), "\n")
  # vector with total numbers of genes controlled by each KD - for normalization
  all_numbers <- lengths(genes_list)
  # collecting KD genes names
  nm_vec <- names(genes_list)
  # creating matrix to populate with numbers of shared genes
  shared_numbers <- matrix(0, nrow = length(genes_list), ncol = length(genes_list))
  rownames(shared_numbers) <- nm_vec
  colnames(shared_numbers) <- nm_vec
  
  
  # populating matrix and list based on directionality and normalization parameters
  for (i in 1:(length(genes_list) - 1)){
    #i <- 2
    cat("\nI'm working with: ", nm_vec[i], "\n", sep = "")
    for (j in (i+1):length(genes_list)){
      #j <- 2
      cat("comparing with ", names(genes_list[j]), "\n")
      colhead <- paste(names(genes_list[i]),"vs", names(genes_list[j]), sep = "")
      shared_only <- genes_list[[i]][genes_list[[i]] %in% genes_list[[j]]]
      all_shared <- length(shared_only)
      cat("Shared only is: ", shared_only, "\n")
      
      
      shared_same_dir <- c()
      shared_diff_dir <- c()
      for(gene in shared_only){
        ratio_i <- ratio_list[[i]][genes_list[[i]] == gene]
        ratio_j <- ratio_list[[j]][genes_list[[j]] == gene]
        print(gene)
        print(ratio_i)
        print(ratio_j)
        
        if(ratio_i > 1 & ratio_j > 1 | ratio_i < 1 & ratio_j < 1){
          # this two genes have same direction
          shared_same_dir <- c(shared_same_dir, gene)
        } else if (ratio_i < 1 & ratio_j > 1 | ratio_i > 1 & ratio_j < 1){
          # this two genes have different direction
          shared_diff_dir <- c(shared_diff_dir, gene)
        }
      }
      shared_same <- length(shared_same_dir)
      shared_diff <- length(shared_diff_dir)
      
      if(all_shared == 0){
        shared_numbers[i,j] <- 0
      } else(
        if(direction == "same"){
          shared_numbers[i,j] <- shared_same
          #cat("shared number is:", shared_numbers[i,j], "\n")
          #shared_numbers[i,j] <- length(shared_same_dir)/mean(c(all_numbers[i], all_numbers[j]))
        } else if(direction == "diff"){
          shared_numbers[i,j] <- shared_diff
        } else if(direction == "all"){
          shared_numbers[i,j] <- all_shared
        }else if (direction == "union"){
          shared_numbers[i,j] <- all_numbers[i] + all_numbers[j] - all_shared
        } else(print("You must specify directionality!"))
      )
    } 
  }
  
  return(shared_numbers)
}

# this function takes ratio matrix and extracts all ratios for all unique 
# gene KD pairs into a vector.
# name_order = forward/reverse - names of combinations can be bind in forward 
# or reversed order (example: gene1vsgene2 or gene2vsgene1)
vectorize_matrix <- function(table, name_order = "forward"){
  all_ratios_vec <- vector(mode = "numeric", length = sum(nrow(table) - seq(1,(ncol(table)-1), 1)))
  all_ratios_names <- vector(mode = "numeric", length = sum(nrow(table) - seq(1,(ncol(table)-1), 1)))
  
  vec_count <- 1
  for(i in 1:(nrow(table)-1)){
    
    for(j in (i+1):(ncol(table))){
      all_ratios_vec[vec_count] <- table[i,j]
      name_i <- colnames(table)[i]
      name_j <- colnames(table)[j]
      if(name_order == "forward"){
        new_name <- paste(name_i, name_j, sep = "vs")
      } else if(name_order == "reverse"){
        new_name <- paste(name_j, name_i, sep = "vs")}
      all_ratios_names[vec_count] <- new_name
      vec_count <- vec_count+1
    }
  }
  names(all_ratios_vec) <- all_ratios_names
  return(all_ratios_vec)
}

# Function: takes a matrix of metrix from normaliz_matrix and transforms
# it into a dataframe with all possible name combinations (therefore
# 36 metrixes become 72 enteries in the dataframe). returns dataframe

kd_metrix_MtoDF <- function(metrics_matrix){
  df_f <- as.data.frame(vectorize_matrix(metrics_matrix, name_order = "forward"))
  colnames(df_f) <- "metrics"
  df_f$comb_names <- row.names(df_f)
  
  df_r <- as.data.frame(vectorize_matrix(metrics_matrix, name_order = "reverse"))
  colnames(df_r) <- "metrics"
  df_r$comb_names <- row.names(df_r)
  
  df <- smartbind(df_f, df_r)
  return(df)
}
#======================================================================
# ALL genes data analysis
#========================================================================

# importing all tables as a list (analysis for all genes detected by RNASeq)
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

# getting list of matrixes generated with different p-value thresholds
matrix_list <- list()
i <- 1
for (p in seq(0.0005, 0.05, by = 0.0005)){
  matrix_list[[i]] <- get_matrix(all_tables_list = all_tables, 
                                                pval_threshold = p, fdr_threshold = F, 
                                                shared_direct = "same", normalized = T)
  i <- i + 1
}


# extracting ratios for each comparison and plotting it in the file:
# creating data frame to populate
kd_names <- row.names(matrix_list[[1]])
comparis_names <- vector()
for(n in kd_names){
  for(n2 in kd_names){
    comparis_names <- append(comparis_names, paste(n, "vs", n2, sep = ""))
  }
}
df_ratios <- data.frame(matrix(ncol = 100, nrow = 81))
row.names(df_ratios) <- comparis_names
colnames(df_ratios) <- seq(0.0005, 0.05, by = 0.0005)

# populating data frame with ratios
count <- 1
for(matrix_ in matrix_list){
  df_ratios[, count] <- as.vector(t(matrix_))
  count <- count +1
}

# building plots and exporting them to the file
pdf(file='ratio_plots.pdf')
for(r in 1:length(row.names(df_ratios))){
  if (sum(df_ratios[r,]) != 0){
    cat("Working with: ",row.names(df_ratios[r,]), "\n")
    plot(x = colnames(df_ratios), y = df_ratios[r, ], xlab = "p-value threshold",
         ylab = "ratio (overlapped/mean_total#)", main = row.names(df_ratios[r,]))
  }
}
dev.off()
###############################################################################

#======================================================================
# Nature Communication genes data anaysis
#========================================================================
# importing all tables as a list (analysis for genes detected by RNASeq AND present
# in Nature communication signature

input_folder <- "./test_files"
all_tables_test <- import_df_asList(input_folder)

test_out <- get_matrix(all_tables_list = all_tables_test, puc = T,
                               pval_threshold = 0.01, fdr_threshold = F, 
                               shared_direct = "all", normalized = F, export = T, 
                               annotation_file_path = "./human_annotation_GTF.txt")

#======================================================================
# Getting ranges of ratios (on demand)
#========================================================================
# getting range of ratios for different parameters (pval, fdr etc.)
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)
test_out <- get_matrix(all_tables_list = all_tables,
                       pval_threshold = 0.05, fdr_threshold = F, 
                       shared_direct = "same", normalized = T, export = T, 
                       annotation_file_path = "./human_annotation_GTF.txt")
max(as.vector(test_out))
all_val <- unique(as.vector(test_out))
all_val <- all_val[-1]
range(all_val)

#=======================================================================
# different normalizations (distribution histograms!)
#========================================================================
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

# plotting as 9 graphs on one page
attach(mtcars)
par(mfrow=c(3,5))
par(mar=c(2,2,2,2))


# getting matrixes with ratios (ATTENTION: you need manualy put p-value 
# threshold in all three normalizations due to loop taking too long to 
# calculate what is needed, and due to the fact that we are doing it just
# for 3 p-val thresholds - so no need to inprove the code to work faster)

pval <- 0.03

testA <- normaliz_matrix(all_tables_list=all_tables, puc = 
                             NULL, pval_threshold = pval, 
                          fdr_threshold = F, normalization = "A")
  
  
testB <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                           pval_threshold = pval, 
                           fdr_threshold = F, normalization = "B")
  
  
testC <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                           pval_threshold = pval, 
                           fdr_threshold = F, normalization = "C")

testD <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "D")

testE <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "E")
  
# extracting values in vector and plotting distribution of values
hist(vectorize_matrix(testA), breaks = 5, ylab = "frequency", xlab = "ratio A",
       main = "")
hist(vectorize_matrix(testB), breaks = 5, ylab = "frequency", xlab = "ratio B",
       main = "")
hist(vectorize_matrix(testD), breaks = 5, ylab = "frequency", xlab = "ratio D",
     main = "")
hist(vectorize_matrix(testE), breaks = 5, ylab = "frequency", xlab = "ratio E",
     main = "")
hist(vectorize_matrix(testC), breaks = 5, ylab = "frequency", xlab = "ratio C",
       main = "")

#========================================================================
# FINDING RESULTS FOR EACH DRIVER IN EACH KD
#========================================================================
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

kd_gene_symbols <- names(all_tables)

# extracting data and exporting tables
new_tables_list <- list()
for (table in 1:length(all_tables)){
  table_name <- names(all_tables[table])
  cat("Wroking with: ", table_name, "\n")
  file_name <- paste(table_name, "_onlyKDgenes.csv", sep = "")
  condition <- all_tables[[table]][["Symbol"]] %in% kd_gene_symbols
  new_table <- all_tables[[table]][eval(condition),]
  write.csv(new_table, file_name, row.names = F)
  new_tables_list[[table_name]] <- new_table
}


# extracting p-values as a matrix. On y-axis (rows) you will see genes that are 
# knocked down, on x-axis (columns) you will see the genes for which we have
# fold change in these KD and p-values for fold change (all these genes are
# the drivers that we used as targets in this project)

first_name <- names(new_tables_list[1])
table_pval <- new_tables_list[[1]][order(new_tables_list[[1]]["Symbol"]), c(1), drop = F]
row.names(table_pval) <- sort(new_tables_list[[1]][,c(7)], decreasing = F)

for(i in 2:length(kd_gene_symbols)){
  first_name[i] <- names(new_tables_list[i])
  table_pval <- cbind(table_pval, 
                      new_tables_list[[i]][order(new_tables_list[[i]]["Symbol"]),+
                                             c(1), drop = F])
}
colnames(table_pval) <- first_name

write.csv(table_pval, "KDdriversInKDexpPval.csv")
# same matrix, but now - fold changes

first_name <- names(new_tables_list[1])
table_FC <- new_tables_list[[1]][order(new_tables_list[[1]]["Symbol"]), c(5), drop = F]
row.names(table_FC) <- sort(new_tables_list[[1]][,c(7)], decreasing = F)

for(i in 2:length(kd_gene_symbols)){
  first_name[i] <- names(new_tables_list[i])
  table_FC <- cbind(table_FC, 
                      new_tables_list[[i]][order(new_tables_list[[i]]["Symbol"]),+
                                             c(5), drop = F])
}
colnames(table_FC) <- first_name
write.csv(table_FC, "KDdriversInKDexpFC.csv")
#===========================================================
# Analysis: proliferation vs # of up genes, and prolif vs # of down genes
#========================================================================
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

up_dn_matrix <- matrix(0, nrow = 2, ncol = length(names(all_tables)))
colnames(up_dn_matrix) <- names(all_tables)
rownames(up_dn_matrix) <- c("up-genes", "dn-genes")
fdr_thr <- 0.1

for(i in 1:length(names(all_tables))){
  up_genes <- sum(all_tables[[i]][["FDR"]] < fdr_thr & all_tables[[i]][["KD/C-all"]] > 1)
  dn_genes <- sum(all_tables[[i]][["FDR"]] < fdr_thr & all_tables[[i]][["KD/C-all"]] < 1)
  
  up_dn_matrix[1,i] <- up_genes
  up_dn_matrix[2,i] <- dn_genes
}

ave_CI_percent <- c("ANP32E"=82.5, "CDCA8"= 63.6, "DTL"= 9.9, 
                    "EXO1"=56.4, "ITGB3BP"=72.1, "KPNA2"=75.2, 
                    "NEK2"=69.3, "S100PBP"=20.0, "TPX2"= 45.4)

plot(x = ave_CI_percent, y=up_dn_matrix[1,])
plot(x = ave_CI_percent, y=up_dn_matrix[2,])
#=========================================================================
# Appending different ratio (metrix) to table of pooled results of how
# each KD regulates other kd-gene-targets (example: we KD NEK2 but want
# to see shat happens with the rest of 8 genes that we use to kd in other
# experiments)

kd_brb_extract <- read.csv("KDbehavInKD.csv", header = T, stringsAsFactors = F, 
                           check.names = F)

input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

pval <- 0.03

testA <- normaliz_matrix(all_tables_list=all_tables, puc = 
                           NULL, pval_threshold = pval, 
                         fdr_threshold = F, normalization = "A")


testB <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "B")


testC <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "C")

testD <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "D")

testE <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "E")


  
test <- kd_metrix_MtoDF(testA)
#===============================================================================
# Collecting all ratios and numbers under different p-value thresholds for 
# all driver combinations:tables, matrix
#===============================================================================
# function for merging ratios into table
m_to_df <-   function(test, name){
  df <- as.data.frame(vectorize_matrix(test))
  df$comp_names <- row.names(df)
  newdf <- df[,c(2,1)]
  colnames(newdf) <- c("comp_names", name)
  return(newdf)
}


input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)


pval <- 0.001

testA <- normaliz_matrix(all_tables_list=all_tables, puc = 
                           NULL, pval_threshold = pval, 
                         fdr_threshold = F, normalization = "A")


testB <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "B")


testC <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "C")

testD <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "D")

testE <- normaliz_matrix(all_tables_list=all_tables, puc = NULL, 
                         pval_threshold = pval, 
                         fdr_threshold = F, normalization = "E")

write.csv(testA, "shared_same.all_shared.csv")
write.csv(testB, "shared_same.union.csv")
write.csv(testC, "shared_same-shared_different.union.csv")
write.csv(testD, "shared_different.union.csv")
write.csv(testE, "all_shared.union.csv")

testAvec <- as.data.frame(vectorize_matrix(testA))
testAvec$comp_names <- row.names(testAvec)
newAll <- testAvec[,c(2,1)]
colnames(newAll) <- c("comp_names", "shared_same/all_shared")
newAll <- merge(newAll, m_to_df(testB, "shared_same/union"))
newAll <- merge(newAll, m_to_df(testC, "(shared_same - shared_different)/union"))
newAll <- merge(newAll, m_to_df(testD, "shared_different/union"))
newAll <- merge(newAll, m_to_df(testE, "all_shared/union"))


same5perc <- shared_numbers(all_tables_list=all_tables, puc = NULL, 
               pval_threshold = pval, fdr_threshold = F, direction = "same")
newAll <- merge(newAll, m_to_df(same5perc, "same_dir#"))

diff5perc <- shared_numbers(all_tables_list=all_tables, puc = NULL, 
                            pval_threshold = pval, fdr_threshold = F, 
                            direction = "diff")
newAll <- merge(newAll, m_to_df(diff5perc, "diff_dir#"))

all5perc <- shared_numbers(all_tables_list=all_tables, puc = NULL, 
                            pval_threshold = pval, fdr_threshold = F, 
                            direction = "all")
newAll <- merge(newAll, m_to_df(all5perc, "allshared#"))

union5perc <- shared_numbers(all_tables_list=all_tables, puc = NULL, 
                            pval_threshold = pval, fdr_threshold = F, 
                            direction = "union")
newAll <- merge(newAll, m_to_df(union5perc, "union#"))

write.csv(newAll, "allRatios&NumbersTable.csv", row.names = F)

#===============================================================================
# PUC for GE vs proliferation and KD GE change - creating table with genes
# that pass puc under each KD
#===============================================================================
#extra function for this task:
exFullsubtab <- function(all_tables, GEvsP, extract){
  if(extract == "fc"){
    a <- 5
  } else if (extract == "pval"){
    a <- 1
  }
  table <- all_tables[[1]][all_tables[[1]][["UniqueID"]] %in% GEvsP$UniqueID, c(6,a)]
  colnames(table) <- c("UniqueID", names(all_tables[1]))
  
  for(i in 2:length(names(all_tables))){
    cat("Processing table: ", names(all_tables[i]), "\n")
    temp_table <- all_tables[[i]][all_tables[[i]][["UniqueID"]] %in% GEvsP$UniqueID, c(6,a)]
    colnames(temp_table) <- c("UniqueID", names(all_tables[i]))
    table <- merge(table, temp_table)
    rm(temp_table)
  }
  return(table)
}
####

input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

# We decided to stick to ge vs proliferation analysis #3: correlation dir
# checked in 3 groups (all, cntr, kd), fisher and fisher_fdr calculated only
# on cntr and kd groups.
# Importing GE vs Proliferation analysis table
GEvsP <- read.csv("123dir23FishFdr.csv", header = T, check.names = F)
GEvsP$GEvsPcor <- sign(GEvsP$all_rho)

fc_table <- exFullsubtab(all_tables = all_tables, GEvsP = GEvsP, extract = "fc")
pv_table <- exFullsubtab(all_tables = all_tables, GEvsP = GEvsP, extract = "pval")
#sum(GEvsP$UniqueID == pv_table$UniqueID)

# extracting pvalues for the genes that pass puc criteria
ext_vecTem <- fc_table[2] > 1 & GEvsP$GEvsPcor == -1 | fc_table[2] < 1 & GEvsP$GEvsPcor == 1
newGEvsPTem <- merge(GEvsP, pv_table[ext_vecTem,c(1,2)], all = T)
rm(ext_vecTem)
for(i in 3:10){
  ext_vecTem <- fc_table[i] > 1 & GEvsP$GEvsPcor == -1 | fc_table[i] < 1 & GEvsP$GEvsPcor == 1
  newGEvsPTem <- merge(newGEvsPTem, pv_table[ext_vecTem,c(1,i)], all = T)
  rm(ext_vecTem)
}

# replacing NA with 1 - if puc is not met, we make pvalue = 1
test <- newGEvsPTem[,19:27]
test[is.na(test)] <- 1
newGEvsPTem[,19:27] <- test
colnames(newGEvsPTem)[19:27] <- as.vector(paste(names(newGEvsPTem[,19:27]), "pval", sep = "_"))
  
# exporting table
write.csv(newGEvsPTem, "GEvsP_PUCpass.csv", row.names = F)
#===============================================================================
# Estimate ratio range for each KD-target
#===============================================================================
source("./cc_drivers_functions/ratioStat.R")
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)


matrixO1 <- normaliz_matrix(all_tables_list=all_tables, puc = 
                           NULL, pval_threshold = 0.01, 
                         fdr_threshold = F, normalization = "B")
matrixOO1 <- normaliz_matrix(all_tables_list=all_tables, puc = 
                               NULL, pval_threshold = 0.001, 
                             fdr_threshold = F, normalization = "B")


statO1 <- ratioStat(matrixO1, 0.01)
statOO1 <- ratioStat(matrixOO1, 0.001)
#===============================================================================
# TEMP
#===============================================================================
union5perc <- shared_numbers(all_tables_list=all_tables, puc = NULL, 
                             pval_threshold = 0.03, fdr_threshold = F, 
                             direction = "union")

unionvecF <- vectorize_matrix(union5perc, name_order = "forward")
unionvecR <- vectorize_matrix(union5perc, name_order = "reverse")
allunion <- c(unionvecF, unionvecR)
write.csv(allunion,"temp.csv")

#####
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)
pval_table <- all_tables[[1]][,c(6,1)]
colnames(pval_table) <- c("UniqueID", names(all_tables)[1])

for(i in 2:length(all_tables)){
  cat("Processing: ", names(all_tables)[i], "\n")
  new_table <- all_tables[[i]][,c(6,1)]
  new_t_names <- c("UniqueID", names(all_tables)[i])
  colnames(new_table) <- new_t_names
  pval_table <- merge(pval_table, new_table, all = T)
}

pval_table$min_pval <- apply(pval_table[, c(2:10)], 1, function(x) min(x))

write.csv(pval_table, "brb-pval-forKDtable.csv", row.names = F)
