# This is the script for outputting a rudimetary KO Hierarchy & Counts Tables

# Install and load packages
install.packages("jsonlite")
install.packages("stringr")
install.packages("dplyr")
install.packages("readr")
library(jsonlite)
library(stringr)
library(dplyr)
library(readr)

# Custom Function definitions

# Function "create_KO_mapfile" creates a dataframe version of the KEGG Ontology Database
# from the provided JSON file.
# For each gene node in the JSON file, the function finds the path (in terms of levels)
# to get to said node. It then appends said path to a dataframe 
# Parameters: json_file = KO JSON database file

create_KO_mapfile <- function(json_file) {
  # Require jsonlite package for the function
  require(jsonlite)
  
  # Read the json file from the provided filename in the parameters
  node <- read_json(json_file)
  
  # Initialize mapfile dataframe, path vector.
  mapfile <- data.frame()
  path <- c()
  
  # This loop goes through every nested category in the json file
  function_index <- 1
  
  while (function_index < (length(node$children) + 1)){
    subfunction_index <- 1
    
    while (subfunction_index < (length(node$children[[function_index]]$children) + 1)){
      pathway_index <- 1
      
      while (pathway_index < (length(node$children[[function_index]]
                                     $children[[subfunction_index]]$children) + 1)){
        gene_index <- 1
        
        while (gene_index < (length(node$children[[function_index]]
                                    $children[[subfunction_index]]
                                    $children[[pathway_index]]$children) + 1)){
          
          # This part of the loop retrieves the path required to reach every individual gene entry in
          # the JSON file. The path consists of four categories/levels, from function to subfunction
          # to pathway to gene.
          
          path <- c(str_sub(node$children[[function_index]]$name, 7),
                    str_sub(node$children[[function_index]]$children[[subfunction_index]]$name, 7),
                    str_sub(node$children[[function_index]]$children[[subfunction_index]]$children[[pathway_index]]$name, 7),
                    str_sub(node$children[[function_index]]$children[[subfunction_index]]
                            $children[[pathway_index]]$children[[gene_index]]$name, 9),
                    str_sub(node$children[[function_index]]$children[[subfunction_index]]
                            $children[[pathway_index]]$children[[gene_index]]$name, end = 6))
          
          # FOR TROUBLESHOOTING: Print the current value of the indices
          print(paste(function_index, subfunction_index, pathway_index, gene_index))
          
          # Append the path vector to the mapfile dataframe
          mapfile <- rbind(mapfile, path)
          
          gene_index <- gene_index + 1
        }
        pathway_index <- pathway_index + 1
        
      }
      subfunction_index <- subfunction_index + 1
    }
    
    function_index <- function_index + 1
  }
  
  return(mapfile)
}

# The function "summarize_KO" produces a data frame that summarizes the expression levels
# of each function/subfunction/pathway/gene. This assumes that the KO codes are the row names!
# This also uses the previous function.
# Parameters: mapfile = KEGG Orthology mapfile;
#            counts = counts from KO_metagenome_out/pred_metagenome_unstrat.tsv
#            level = hierarchy level/category to generate. Must be an integer from 1 to 4
#                     Default is 3 (L3; pathway)
#            limit = number of KO codes/genes to process. Only processes the first n valid
#                     (i.e. present in the mapfile) KO codes in the counts table
#                     More useful for troubleshooting. NULL by default (no limit)
# Outputs: a data frame. You may save the data frame to an object by
# assigning the function call to said object (e.g. df <- summarize_KO())

summarize_KO <- function(mapfile, counts, level = 3, limit = NULL){
  # Require dplyr for this function (function uses cbind, rbind, pipeline operator)
  require(dplyr)
  
  # Initialize index (to be used if parameter 'limit' is defined) and dataframe
  # for containing the generated hierarchy table
  i <- 1
  df <- data.frame()
  
  # This loop goes through every row name (i.e. KO code) in the counts table
  for (row in rownames(counts)){
    
    # FOR TROUBLESHOOTING: Print row (may also serve as rudimentary progress indicator)
    # Remove '#' if to be used
    # print(row)
    
    # This nested if/else determines if the current KO code (in variable "row") is present in the
    # mapfile. If present in the mapfile, it will proceed to retrieve the categories in which said KO
    # code is present. If absent from the mapfile, it will notify the user through a warning,
    if (any(sapply(mapfile$KO, function(x){x == row}, simplify = TRUE, USE.NAMES = FALSE))){
      if (level == 1){ # For Level 1: function
        l1_list<- data.frame(Function = mapfile[mapfile$KO == row, 1])
        l1_counts <-  counts[rep(which(rownames(counts) == row), nrow(l1_list)), ]
        l1_list <- cbind(l1_list, l1_counts)
        
        df <- rbind(df, l1_list)
        
        colnames(df)[1] <- "Function"
      }
      else if (level == 2){ # For Level 2: Subfunction
        l2_list <- data.frame(Subfunction = mapfile[mapfile$KO == row, 2])
        l2_counts <-  counts[rep(which(rownames(counts) == row), nrow(l2_list)), ]
        l2_list <- cbind(l2_list, l2_counts)
        
        df <- rbind(df, l2_list)
        
        colnames(df)[1] <- "Subfunction"
      }
      else if (level == 3){ # For Level 3: Pathway
        l3_list <- data.frame(Pathway = mapfile[mapfile$KO == row, 3])
        l3_counts <-  counts[rep(which(rownames(counts) == row), nrow(l3_list)), ]
        l3_list <- cbind(l3_list, l3_counts)
        
        df <- rbind(df, l3_list)
        
        colnames(df)[1] <- "Pathway"
      }
      else if (level == 4){ # For Level 4: Gene. Only adds a description.
        l4_desc <- mapfile[mapfile$KO == row, 4][1]
        df <- rbind(df, cbind(Description = l4_desc, KO = row, counts[rownames(counts) == row,]))
      }
      else{ # If level is incorrectly defined
        stop("You have entered an invalid value for parameter: level. Please provide a number from 1 to 4.")
      }
    }
    else{
      message(paste(row, " does not exist in the KEGG Orthology BRITE Hierachies.
                  Make sure to remove this from your data"))
      i <- i - 1
    }
    
    # If parameter "limit" is defined, this part stops the loop once said limit is reached
    i <- i + 1
    if (!is.null(limit)){
      if(i > limit){
        break
      }
    }
  }
  
  # Sum the counts of each category appended in the dataframe.
  # Skipped if level = 4 (gene level does not need summarizing)
  if(level != 4){
    df <- df %>% group_by_at(1) %>% summarize_all(sum)
  }
  
  return(df)
}

# The function map_KO generates a mapfile from gene to specified KO level
# and is essentially a modified version of the summarize_KO() function, just
# without making a counts summary.
# All parameters are the same as that of summarize_KO

map_KO <- function(mapfile, counts, level = 3, limit = NULL){
  # Require dplyr for this function (function uses cbind and rbind)
  require(dplyr)
  
  # Initialize index (to be used if parameter 'limit' is defined) and dataframe
  # for containing the geenrated hierarchy table
  i <- 1
  df <- data.frame()
  
  # This loop goes through every row name (i.e. KO code) in the counts table
  for (row in rownames(counts)){
    # FOR TROUBLESHOOTING: Print row (may also serve as rudimentary progress indicator)
    # Remove '#' if to be used
    # print(row)
    
    # This nested if/else determines if the current KO code (in variable "row") is present in the
    # mapfile. If present in the mapfile, it will proceed to retrieve the categories in which said KO
    # code is present. If absent from the mapfile, it will notify the user through a warning,
    if (any(sapply(mapfile$KO, function(x){x == row}, simplify = TRUE, USE.NAMES = FALSE))){
      if (level == 1){ # For Level 1: function
        l1_list<- data.frame(Function = mapfile[mapfile$KO == row, 1], Subfunction = mapfile[mapfile$KO == row, 2],
                             Pathway = mapfile[mapfile$KO == row, 3], Gene = mapfile[mapfile$KO == row, 4], KO = mapfile[mapfile$KO == row, 5])
        l1_counts <-  counts[rep(which(rownames(counts) == row), nrow(l1_list)), ]
        l1_list <- cbind(l1_list, l1_counts)
        
        df <- rbind(df, l1_list)
        
        colnames(df)[1] <- "Function"
      }
      else if (level == 2){ # For Level 2: Subfunction
        l2_list <- data.frame(Subfunction = mapfile[mapfile$KO == row, 2],
                              Pathway = mapfile[mapfile$KO == row, 3], Gene = mapfile[mapfile$KO == row, 4], KO = mapfile[mapfile$KO == row, 5])
        l2_counts <-  counts[rep(which(rownames(counts) == row), nrow(l2_list)), ]
        l2_list <- cbind(l2_list, l2_counts)
        
        df <- rbind(df, l2_list)
        
        colnames(df)[1] <- "Subfunction"
      }
      else if (level == 3){ # For Level 3: Pathway
        l3_list <- data.frame(Pathway = mapfile[mapfile$KO == row, 3], Gene = mapfile[mapfile$KO == row, 4], KO = mapfile[mapfile$KO == row, 5])
        l3_counts <-  counts[rep(which(rownames(counts) == row), nrow(l3_list)), ]
        l3_list <- cbind(l3_list, l3_counts)
        
        df <- rbind(df, l3_list)
        
        colnames(df)[1] <- "Pathway"
      }
      else if (level == 4){ # For Level 4: Gene. Only adds a description.
        l4_desc <- mapfile[mapfile$KO == row, 4][1]
        df <- rbind(df, cbind(Description = l4_desc, KO = row, counts[rownames(counts) == row,]))
      }
      else{ # If level is incorrectly defined
        stop("You have entered an invalid value for parameter: level. Please provide a number from 1 to 4.")
      }
    }
    else{
      # This eliminates KO entries that are not present in the mapfile (as they cannot be counted)
      message(paste(row, " does not exist in the KEGG Orthology BRITE Hierachies.
                  Make sure to remove this from your data"))
      i <- i - 1
    }
    
    # If parameter "limit" is defined, this part stops the loop once said limit is reached
    i <- i + 1
    if (!is.null(limit)){
      if(i > limit){
        break
      }
    }
  }
  
  return(df)
}

# The function toLog transforms all numeric data in a data frame to a logarithm
# Default logarithm is base 10; type 'exp(1)' if you want the natural logarithm

toLog <- function(df, base = 10, add_constant){
  row_limit <- length(rownames(df))
  col_limit <- length(colnames(df))
  row_index <- 1
  
  while (row_index < row_limit+1) {
    col_index <- 1
    
    while (col_index < col_limit+1){
      if(is.numeric(df[row_index, col_index])){
        df[row_index, col_index] <- log((df[row_index, col_index] + add_constant), base)
      }
      
      col_index <- col_index + 1
    }
    row_index <- row_index + 1
  }
  
  return(df)
}


# KO categorization begins here

# Read the JSON file for the KEGG BRITE hierarchies and transform into a dataframe. 
# Extract this file from the mapfiles folder
# This will serve as the KO database
# This JSON file is sourced from https://www.genome.jp/kegg-bin/get_htext?query=&htext=ko00001
json_file <- "kegg_britemap.json"
KO_mapfile <- create_KO_mapfile(json_file)

# OPTIONAL: If you have the mapfile in tsv form already (because generating the mapfile from scratch
# takes a while), read it through the following line:
KO_mapfile <- read_tsv("<path_to_tsv_mapfile>")

KO_mapfile <- read_tsv("/Users/regaladoricryan/Downloads/KEGG Brite - 9-30-24.tsv")

# This will prune out genes relating to human diseases and organismal systems that is irrelevant for downstream analyses.
KO_mapfile_sans_human <- subset(KO_mapfile, !((Function == "Human Diseases") | (Function == "Organismal Systems")
                                              |(Function == "Brite Hierarchies") | (Function == "Not Included in Pathway or Brite")))
KO_mapfile_sans_human <- subset(KO_mapfile_sans_human, !(grepl("eukaryote", KO_mapfile_sans_human$Subfunction)|
                                                           grepl("eukaryote", KO_mapfile_sans_human$Pathway)))

# This reads the pred metagenome unstrat file and retrieves the counts of each KO per sample
KO_counts <- read.table("picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv", header=TRUE, row.names=1,
                        check.names=FALSE, sep="\t")
rownames(KO_counts) <- str_sub(rownames(KO_counts), start = 4)

# This creates the hierarchy data from the mapfile without human genes
KO_hierarchy_func <- summarize_KO(KO_mapfile_sans_human, KO_counts, level = 1)
KO_hierarchy_subf <- summarize_KO(KO_mapfile_sans_human, KO_counts, level = 2)
KO_hierarchy_path <- summarize_KO(KO_mapfile_sans_human, KO_counts, level = 3)
KO_hierarchy_desc <- summarize_KO(KO_mapfile_sans_human, KO_counts, level = 4)

# Optional: Remove "Others" category from Level 3 (Pathways)
KO_hierarchy_path <- subset(KO_hierarchy_path, !(Pathway == "Others"))

# Export the individual hierarchies/counts table
write.table(KO_hierarchy_func, "<output_path_func.tsv>", sep = "\t", quote=FALSE, row.names=FALSE)
write.table(KO_hierarchy_subf, "<output_path_subf.tsv>", sep = "\t", quote=FALSE, row.names=FALSE)
write.table(KO_hierarchy_path, "<output_path_path.tsv>", sep = "\t", quote=FALSE, row.names=FALSE)
write.table(KO_hierarchy_desc, "<output_path_gene.tsv>", sep = "\t", quote=FALSE, row.names=FALSE)