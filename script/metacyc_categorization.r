# This is the script for outputting a rudimetary MetaCyc Hierarchy & Counts Table

# Install and load the required package
install.packages("stringr")
install.packages("tibble")
library(stringr)
library(tibble)

# Custom Function definitions

# Function "get_MC_levels" extracts all parents (e.g. higher levels containing a particular MC pathway)
# It then ranks these parents according to their level in the MetaCyc Pathway Hierarchy
# It then assigns a Level I parent and a Level II Parent to the pathway of interest

# PARAMETERS: pathway_desc = pathway description of interest, extracted by create_MC_hierarchies();
#           MC_lv1_children = database file for Level 1 parents and their Level 2 children;
#           MC_path_parents = database file for pathways and their parent/s;
#           MC_file_separator = file separator for the databases. " // " by default;
#           makeDescriptionAsDefaultL2 = if TRUE, assign pathway description also as the L2
#                           if FALSE, assign L1 parent also as L2 parent. This is passed to first_L2()

get_MC_levels <- function(pathway_desc, MC_lvl1_children, MC_path_parents, MC_file_separator, makeDescriptionAsDefaultL2) {
  cur_MC_levels <- data.frame()
  original_pathway_desc <- pathway_desc
  
  # Remove quotation marks
  pathway_desc <- gsub('"', '', pathway_desc)
  
  # Check if the pathway description is blank or NA
  if(((identical(pathway_desc,""))|(is.na(pathway_desc)))|(is.null(pathway_desc))){
    message(paste0("No pathway is provided. Please check your data."))
    return(c(NA, NA))
  }
  else{
    
    # This deals with parentheses in the pathway description
    # This encodes the escape character for the parenthesis, '\\(', into the pathway description
    # which allows for string matching later on
    if(grepl("\\(", pathway_desc) | grepl("\\)", pathway_desc)){
      pathway_desc <- gsub("\\(", "\\\\(", pathway_desc)
      pathway_desc <- gsub("\\)", "\\\\)", pathway_desc)
    }
    
    # This checks if the pathway description matches any of those in MC_path_parents
    if(!identical(grep(pathway_desc, MC_path_parents$Pathways, ignore.case = TRUE), integer(0))){
      
      # If successful, it takes the parents of the matching pathway description
      for(index in as.numeric(grep(pathway_desc, MC_path_parents$Pathways, ignore.case = TRUE))){
        parents1 <- strsplit(MC_path_parents[[2]][index], MC_file_separator)
        parents2 <- strsplit(MC_path_parents[[3]][index], MC_file_separator)
        parents3 <- strsplit(MC_path_parents[[4]][index], MC_file_separator)
        parents4 <- strsplit(MC_path_parents[[5]][index], MC_file_separator)
        parents5 <- pathway_desc
      }
      
      path_parents <- data.frame(Parents = c(unlist(parents1), unlist(parents2), unlist(parents3), unlist(parents4)))
      path_parents <- subset(path_parents, path_parents[[1]] != "Pathways")
      path_parents <- subset(path_parents, path_parents[[1]] != "")
      path_parents <- subset(path_parents, !(grepl("  ", path_parents[[1]])))
      
      # This feeds the list of parents into rank_MC_levels() for ranking
      path_parents_ranked <- rank_MC_levels(path_parents, MC_lvl1_children)
      
      # This feeds the rank 1 (i.e. level 1) parents into first_L1() for L1 assignment
      # If there are no rank 1s listed initially, L1s are identified for each L2 parent present
      path_parents_rank1 <- subset(path_parents_ranked, Rank == 1)
      if (nrow(path_parents_rank1) == 0){
        path_parents_rank1 <- find_L1(subset(path_parents_ranked, Rank == 2), MC_lvl1_children)} 
      first_rank1_parent <- first_L1(path_parents_rank1, MC_lvl1_children)
      
      #This feeds the rank 1 (i.e. level 1) parents into first_L2() for L2 assignment
      path_parents_rank2 <- subset(path_parents_ranked, Rank == 2)
      first_rank2_parent <- first_L2(path_parents_rank2, MC_lvl1_children, first_rank1_parent,
                                     original_pathway_desc, makeDescriptionAsDefaultL2)
      
      return(c(first_rank1_parent, first_rank2_parent))
    }
    else{
      message(paste0("The pathway \"", pathway_desc, "\" does not exist in the MetaCyc Pathway Hierarchies.
                  Make sure to remove this from your data"))
      
      return(c(NA, NA))
    }
  }
}


# function "rank_MC_levels" ranks all parents extracted by "get_MC_levels()"
# Parents are compared to the L1/L2 database, and assigned rank 1 or 2 based on which level they match
# Other parents are assigned L3, and are deprioritized.
# PARAMETERS: path_parents = pathway parents' dataframe;
#             MC_lv1_children = database file for Level 1 parents and their Level 2 children

rank_MC_levels <- function(path_parents, MC_lvl1_children){
  ranks <- data.frame()
  
  for(parent in path_parents[[1]]){
    isL1 <- any(grepl(parent, MC_lvl1_children$Pathways))
    isL2 <- any(grepl(parent, MC_lvl1_children$Subclasses))
    isL3 <- !(isL1 | isL2)
    
    if(isL1 == TRUE){
      parent_rank <- 1
    }
    else if(isL2 == TRUE){
      parent_rank <- 2
    }
    else if(isL3 == TRUE){
      parent_rank <- 3
    }
    else{
      parent_rank <- NA
    }
    
    ranks <- rbind(ranks, parent_rank)
  }
  
  path_parents <- cbind(path_parents, ranks)
  names(path_parents)[2] <- "Rank"
  return(path_parents)
}

# function first_L1 and first_L2 retrieve the first Level 1 and Level 2 parent that appears on the
# L1/L2 database. MetaCyc pathways may have multiple L1 or L2 parents, and this resolves the conflict
# by assigning it the "more important" L1/L2 parent

# PARAMETERS: path_parents_rank1 = rank 1 pathway parents' dataframe;
#           path_parents_rank2 = rank 2 pathway parents' dataframe;
#           MC_lv1_children = database file for Level 1 parents and their Level 2 children
#           L1 parent = L1 parent assigned to the pathway. This is used to get Level 2.
#           sortOthersLast = sorts "Others..." categories last, even if they are ranked higher in
#           the metacyc hierarchies. "Others" will only be assigned if a pathway has no other
#           L2 candidates. For first_L2() only.
#           makeDescriptionAsDefaultL2 = if TRUE, assign pathway description also as the L2
#                           if FALSE, assign L1 parent also as L2 parent. For first_L2 only

first_L1 <- function(path_parents_rank1, MC_lvl1_children){
  cur_first_L1 <- NULL
  cur_first_L1_index <- NULL
  
  for(parent in path_parents_rank1[[1]]){
    if (any(grepl(parent, MC_lvl1_children$Pathways))){
      parent_index <- as.numeric(grep(parent, MC_lvl1_children$Pathways, fixed = TRUE))
      
      for (candidate_index in parent_index){
        if (MC_lvl1_children$Pathways[candidate_index] == parent){
          parent_index <- candidate_index
        }
      }
      
      
    }
    else{
      parent_index <- NA
    }
    
    if (is.na(parent_index)){
      next
    }
    else if (is.null(cur_first_L1)){
      cur_first_L1 <- parent
      cur_first_L1_index <- parent_index
    }
    else if (parent_index < cur_first_L1_index){
      cur_first_L1 <- parent
      cur_first_L1_index <- parent_index
    }
  }
  
  return(cur_first_L1)
}

first_L2 <- function(path_parents_rank2, MC_lvl1_children, L1_parent, pathway_desc,
                     sortOthersLast = TRUE, makeDescriptionAsDefaultL2 = TRUE){
  cur_first_L2 <- NULL
  cur_first_L2_index <- NULL
  other_L2 <- ""
  cur_first_L2_index <- NULL
  
  if(identical(MC_lvl1_children$Subclasses[MC_lvl1_children$Pathways == L1_parent],"")){
    if(makeDescriptionAsDefaultL2 == TRUE){
      cur_first_L2 <- paste(pathway_desc, "(L2)")
    }
    else{
      cur_first_L2 <- paste(L1_parent, "(L2)")
    }
  }
  else{
    for(parent in path_parents_rank2[[1]]){
      
      if (grepl(gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", parent)),
                MC_lvl1_children$Subclasses[MC_lvl1_children$Pathways==L1_parent])){
        parent_index <- as.numeric(regexpr(parent,
                                           MC_lvl1_children$Subclasses[MC_lvl1_children$Pathways==L1_parent], fixed = TRUE))
      }
      else{
        parent_index <- NA
      }
      
      if (is.na(parent_index)){
        next
      }
      else {
        if (is.null(cur_first_L2)){
          cur_first_L2 <- parent
          cur_first_L2_index <- parent_index
        }
        else{
          if(sortOthersLast == TRUE){
            if(grepl("other", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", parent)), ignore.case = TRUE)){
              other_L2 <- parent
              other_L2_index <- parent_index
            }
            else if (parent_index < cur_first_L2_index){
              cur_first_L2 <- parent
              cur_first_L2_index <- parent_index
            }
          }
          else{
            if (parent_index < cur_first_L2_index){
              cur_first_L2 <- parent
              cur_first_L2_index <- parent_index
            }
          }
        }
        
      }
    }
    
    if(sortOthersLast == TRUE){
      if((is.null(cur_first_L2))&(grepl("other", other_L2, ignore.case = TRUE))){
        cur_first_L2 <- other_L2
      }
    }
    
    # If the function somehow goes through all L2 children but still find nothing, this assigns
    # either pathway name or L1 parent as the L2 level
    if (is.null(cur_first_L2)){
      if(makeDescriptionAsDefaultL2 == TRUE){
        cur_first_L2 <- paste(pathway_desc, "(L2)")
      }
      else{
        cur_first_L2 <- paste(L1_parent, "(L2)")
      }
    }
  }
  
  
  
  return(cur_first_L2)
}

# function find_L1 retrieves the Level 1 parents associated with any Level 2 parent associated
# with a given pathway. Note that this function will only be used in the case that there are no L1
# parents in the pathway and parents database (i.e. MC_path_parents); however, L2 parents MUST
# be present. If no L1 or L2 parents are identified, adjust the mapfile for MC_path_parents.

# PARAMETERS: path_parents_rank2 = rank 2 pathway parents' dataframe;
#           MC_lv1_children = database file for Level 1 parents and their Level 2 children

find_L1 <- function(path_parents_rank2, MC_lvl1_children){
  found_L1_parents <- data.frame()
  
  for(L1_parent in MC_lvl1_children$Pathways){
    for(parent in path_parents_rank2[[1]]){
      if (grepl(gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", parent)),
                MC_lvl1_children$Subclasses[MC_lvl1_children$Pathways==L1_parent])){
        found_L1_parents <- rbind(found_L1_parents, c(L1_parent, 1))
      }
    }
  }
  
  names(found_L1_parents)[1] <- "Parents"
  names(found_L1_parents)[2] <- "Ranks"
  
  return(found_L1_parents)
}

# The function "create_MC_hierarchy" produces a data frame that summarizes the hierarchy
# of each MetaCyc pathway.This assumes that the MC pathway descriptions are present!

# Parameters: node = hierarchy; counts = counts from KO_metagenome_out/pred_metagenome_unstrat.tsv

# PARAMETERS:MC_lv1_children = database file for Level 1 parents and their Level 2 children;
#           MC_path_parents = database file for pathways and their parent/s;
#           MC_abundance_w_desc. Pathway abundance file ("path_abun_unstrat.tsv.gz by default)
#           MC_file_separator = file separator for the databases. " // " by default;
#           makeDescriptionAsDefaultL2 = if TRUE, assign pathway description also as the L2
#                           if FALSE, assign L1 parent also as L2 parent. This is passed to first_L2()
#           runLimit = if assigned a number, create MC hierarchies only for the specified number of pathways
#                           default is NULL; MC hierarchies will be created for all pathways listed
#                           Useful for diagnosis

# Outputs: a data frame. You may save the data frame to an object by
# assigning the function call to said object (e.g. df <- create_KO_hierarchy())

create_MC_hierarchy <- function(MC_lvl1_children, MC_path_parents, MC_abundance_w_desc,
                                MC_file_separator = " // ", makeDescriptionAsDefaultL2 = TRUE, runLimit = NULL){
  i <- 1
  df <- data.frame()
  
  for (desc in MC_abundance_w_desc[[2]]){
    path <- get_MC_levels(desc, MC_lvl1_children, MC_path_parents, MC_file_separator, makeDescriptionAsDefaultL2)
    
    if (identical(path,c(NA,NA))){
      for(j in 1:2){
        df[i, j] <- NA
      }
    }
    else{
      for(j in 1:2){
        df[i, j] <- path[j]
      }
    }
    
    i <- i + 1
    
    if (!is.null(runLimit)){
      if (i > runLimit) { return(df) }
    }
  }
  
  return(df)
}

# The function "add_MC_descriptions" adds the pathway descriptions to the counts data table. This is
# only reuqired if the descriptions were not generated directly from PICRUST2

add_MC_descriptions <- function(MC_pathways, MC_descriptions, abundanceFileInput = FALSE){
  description_df <- data.frame()
  i <- 1
  
  if (abundanceFileInput == TRUE){
    for(extpath in rownames(MC_abundance)){
      if(!is.na(extpath)){
        index <- grep(sub("\\(", "\\\\(", sub("\\)", "\\\\)", extpath)), MC_descriptions$Pathway, fixed = TRUE)
        
        if(!identical(index, integer(0))){
          description_df <- rbind(description_df, MC_descriptions$'Description'[as.numeric(index)])
        }
        else{
          description_df <- rbind(description_df, NA)
        }
      }
      else{
        description_df <- rbind(description_df, NA)
      }
      
      message(paste("Process ongoing.", i, "out of", length(rownames(MC_abundance)), "completed"))
      i <- i + 1
    }
    
    merged_MC_pathways <- cbind(rownames(MC_pathways), description_df, MC_pathways)
    names(merged_MC_pathways)[1] <- "MC Pathway"
    names(merged_MC_pathways)[2] <- "MC Pathway Description"
    return(merged_MC_pathways)
  }
  else{
    for(extpath in MC_pathways$`MC Pathways`){
      if(!is.na(extpath)){
        index <- grep(extpath, MC_descriptions$`Pathway`, fixed = TRUE)
        
        if(!identical(index, integer(0))){
          description_df <- rbind(description_df, MC_descriptions$'Description'[as.numeric(index)])
        }
        else{
          description_df <- rbind(description_df, NA)
        }
      }
      else{
        description_df <- rbind(description_df, NA)
      }
      
      message(paste("Process ongoing.", i, "out of", length(MC_pathways$`MC Pathways`), "completed"))
      i <- i + 1
    }
    
    merged_MC_pathways <- cbind(MC_pathways, description_df)
    colnames(merged_MC_pathways) <- c("MC Pathway", "MC Pathway Description")
    return(merged_MC_pathways)
  }
  
}

# The function "fuseMCPathAndDesc" fuses the pathway codes to the pathway descriptions
# This function, by default, assumes that the pathway codes are provided first, then the descriptions second

fuseMCPathAndDesc <- function(MC_hierarchy_no_NA, pathwayCodeFirst = TRUE){
  L1_and_L2 <- MC_hierarchy_no_NA[,1:2]
  count_data <- MC_hierarchy_no_NA[,5:length(MC_hierarchy_no_NA)]
  
  pathWithDesc <- data.frame()
  
  index_limit <- 1 + length(rownames(MC_hierarchy_no_NA))
  i <- 1
  
  if (pathwayCodeFirst == TRUE){
    while (i < index_limit){
      fusedPathAndDesc <- paste0(MC_hierarchy_no_NA[i,3], ": ", MC_hierarchy_no_NA[i,4])
      pathWithDesc <- rbind(pathWithDesc, fusedPathAndDesc)
      
      i <- i + 1
    }
  }
  else if (pathwayCodeFirst == FALSE){
    while (i < index_limit){
      fusedPathAndDesc <- paste0(MC_hierarchy_no_NA[i,4], ": ", MC_hierarchy_no_NA[i,3])
      pathWithDesc <- rbind(pathWithDesc, fusedPathAndDesc)
      
      i <- i + 1
    }
  }
  else{
    stop("You have specified an invalid input for parameter: pathwayCodeFirst\nParameter may only be TRUE or FALSE")
  }
  
  modified_df <- cbind(L1_and_L2, pathWithDesc, count_data)
  names(modified_df)[3] <- 'MC Pathway'
  return(modified_df)
}

# The function "chooseTopN" chooses the top N pathways from the provided MC hierarchy + counts dataframe
# The function may use a specified algorithm [to follow]. Default sorts by the number of reads

chooseTopN <- function(df, n = 25, sortBy = "number of reads", breakTies = TRUE, addParameters = NULL){
  n <- round(n)
  row_index_limit <- length(rownames(df))
  
  first_data_col <- 1
  i <- 1
  while (i < (length(df)+1)){
    if (is.numeric(df[1, i])){
      first_data_col <- i
      break
    }
    i <- i+1
  }
  last_data_col <- length(df)
  
  if(is.numeric(n) & ((n < row_index_limit) & (n > 0))){
    num_reads <- data.frame()
    ranks <- data.frame()
    
    #identify sums
    row_index <- 1
    while(row_index < (row_index_limit + 1)){
      row_sum <- sum(df[row_index, first_data_col:last_data_col])
      cur_row <- data.frame(index = row_index, sum = row_sum)
      num_reads <- rbind(num_reads, cur_row)
      
      row_index <- row_index + 1
    }
    
    ##rank items
    for(cur_sum in num_reads$sum){
      lower_than <- 1
      
      for(other_sum in num_reads$sum){
        if(cur_sum < other_sum){
          lower_than <- lower_than + 1
        }
      }
      
      ranks <- rbind(ranks, rank = lower_than)
    }
    
    if (breakTies == TRUE){
      for(rank_i in 1:length(rownames(ranks))){
        tied_rank_indices <- which(ranks[,1] == rank_i)
        
        if (length(tied_rank_indices)>1){
          tied_rank_i <- 1
          
          while(tied_rank_i < (length(tied_rank_indices)+1)){
            ranks[tied_rank_indices[tied_rank_i],1] <- rank_i + tied_rank_i - 1
            tied_rank_i <- tied_rank_i + 1
          }
        }
      }
    }
    
    #choose top N
    top_pathway_indices <- which((ranks[,1]>=1)&(ranks[,1]<=n))
    top_pathways <- df[top_pathway_indices,]
    
    return(top_pathways)
  }
  else if(!is.numeric(n)){
    stop("Your input for parameter: sortBy is invalid. Please check again.")
  }
  else if(n >= length(rownames(df))) {
    stop("Your input for parameter: n exceeds or equals the length of your dataset. Please check again.")
  }
  else if(!(n > 0)){
    stop("Your input for parameter: n is invalid (cannot be zero or negative). Please check again.")
  }
  else{
    stop("Your input for parameter: sortBy is invalid. Please check again.")
  }
}



# MC categorization begins here
# Read the input pathway abundance and MetaCyc hierarchy files.

# Read the TSV for the MC Level I to Children database (children = lower levels)
MC_lvl1_children <- read.table("metacyc mapfiles/metacyc_pathways_children.txt", header=T,
                               check.names=F, sep="\t", quote = "", fill = TRUE) # Extract the file from the mapfiles folder

# Read the TSV for the MC Pathways to Parents database (parents = higher levels)
MC_path_parents <- read.table("metacyc mapfiles/metacyc_pathways_parents.txt", header=T,
                              check.names=F, sep="\t", quote = "", fill = TRUE) # Extract the file from the mapfiles folder

# This reads the pathway abundance file and retrieves the counts of each MC pathway per sample
# Do ensure that the second column contains the descriptions
MC_abundance_w_desc <- read.table("picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv", header=T,
                                  check.names=F, sep="\t", quote = "") # removed row.names=1, added quote = ""
names(MC_abundance_w_desc)[1] <- "MC Pathway" 
names(MC_abundance_w_desc)[2] <- "MC Pathway Description" 

# Making the MC Hierarchy

# This creates the hierarchy data frame using the create_MC_hierarchy function.
# Supply separator data for the MC hierarchy files if they are custom. Default is " // "
# Warning: This process is slow
MC_hierarchy <- create_MC_hierarchy(MC_lvl1_children, MC_path_parents, MC_abundance_w_desc)

# This renames the levels of each hierarchy. Feel free to rename these.
colnames(MC_hierarchy) <- c("MetaCyc Level 1", "MetaCyc Level 2")

# This merges the hierarchy data frame with the counts data frame
MC_hierarchy <- cbind(MC_hierarchy, MC_abundance_w_desc)

# Just in case: Store a copy of the MC_hierarchy (though the function finishes quickly anyways)
MC_hierarchy2 <- MC_hierarchy

# Each line removes rows with an instance of NA in the respective columns
# since STAMP will not accept data with NA (STAMP interprets all NAs as the same
# which causes it to interpret data with NAs as not having a strict tree)
MC_hierarchy_no_NA <- subset(MC_hierarchy2, !is.na(`MetaCyc Level 1`))
MC_hierarchy_no_NA <- subset(MC_hierarchy_no_NA, !is.na(`MetaCyc Level 2`))
MC_hierarchy_no_NA <- subset(MC_hierarchy_no_NA, !is.na(`MC Pathway Description`))

# Write the final output files

write.table(MC_hierarchy_no_NA, "<output-hierarchical-pathway-abundance.tsv>", sep = "\t", quote=F, row.names=FALSE)

# If Pathway description and Pathway code need swapping
# Provide the new arrangement of columns in a vector (shown below) IN ORDER
# Swap columns, then export file
new_col_arrangement <- c("MetaCyc Level 1", "MetaCyc Level 2", "MC Pathway Description", "MC Pathway",
                         "sample-name-1a", "sample-name-1b", "sample-name-2a", "sample-name-2b", 
                         "sample-name-3a", "sample-name-3b")

MC_hierarchy_no_NA_swapped <- MC_hierarchy_no_NA[, new_col_arrangement]

# Alternatively, you can use the select() function in the dplyr function
# You can use the selector functions [e.g. contains()]
library(dplyr)
MC_hierarchy_no_NA_swapped <- select(MC_hierarchy_no_NA, c("MetaCyc Level 1","MetaCyc Level 2","MC Pathway Description","MC Pathway",
                                                           contains("<select-colnames>")))

write.table(MC_hierarchy_no_NA_swapped, "<output-hierarchical-pathway-abundance-swapped.tsv>", sep = "\t", quote=F, row.names=FALSE)

# ADDENDUM: Fuse Pathway and Pathway Description

# Fuse Pathway and Pathway Description
MC_hierarchy_fused_Path <- fuseMCPathAndDesc(MC_hierarchy_no_NA)

# Export file
write.table(MC_hierarchy_fused_Path, "<output-hierarchical-pathway-abundance-fused.tsv>", sep = "\t", quote=F, row.names=FALSE)

# Visualize data in heatmap using STAMP software: https://beikolab.cs.dal.ca/software/STAMP