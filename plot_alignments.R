
library(dplyr)
library(Biostrings)
library(pheatmap)
plot_single_subject_alignments<- function(psa, view_indels=FALSE){
  # get insertions or deletions 
  reference_subject <-(strsplit(as.character(pwalign::unaligned(subject(psa))), split = "")[[1]])
  reference_subject <- t(as.matrix(reference_subject))
  
  psa_mat<- as.matrix(psa)
  psa_mat <- rbind(psa_mat, reference_subject)
  
  if(view_indels == TRUE){
    psa_mat<- add_deletions_to_matrix(psa_mat, psa)
    psa_mat <- add_insertions_to_matrix(psa_mat, psa)
  }

  alphabet <- sort(unique(as.vector(psa_mat)))
  alphabet_nums <- seq_len(length(alphabet))  
  alphabet_map <- setNames(alphabet_nums, alphabet)
  psa_mat_num <- matrix(
    as.numeric(alphabet_map[psa_mat]), 
    nrow = nrow(psa_mat),
    ncol = ncol(psa_mat),
    dimnames = dimnames(psa_mat)    
  )
  
  pallette<- "Paired"
  if (length(alphabet)<= RColorBrewer::brewer.pal.info[pallette, "maxcolors"] ){
    base_colors <- RColorBrewer::brewer.pal(length(alphabet), pallette)
  } else {
    base_colors <- rainbow(length(alphabet))
  }
  
  heatmap<-pheatmap(
    psa_mat_num,
    color = base_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend_breaks = alphabet_nums,
    legend_labels = alphabet,
    main = "Alignments"
  )
  return(heatmap)
}


add_insertions_to_matrix<- function(psa_mat, psa){
  insertions<-pwalign::insertion(psa)
  # for each insertion, add a column at that location 
  for(read_index in seq_along(insertions) ){ # insertions for the reads
    read_insertion <- unlist(insertions[read_index,])
    for(insertion_index in seq_along(read_insertion)){ # for each start location for this read
      new_column<- matrix(nrow = nrow(psa_mat), ncol = width(read_insertion[insertion_index,])) # rows are reads
      new_column[read_index, 1:width(read_insertion[insertion_index,]) ] <- "Ins" # record there's an insertion for read number i 
      psa_mat <- cbind(psa_mat[, 1:insertion_index],# place new column at insertion position
                         new_column, 
                         psa_mat[, (insertion_index + 1):ncol(psa_mat)])
    }
  }
  psa_mat[is.na(psa_mat)] <- "None"
  return(psa_mat)
}

add_deletions_to_matrix<- function(psa_mat, psa){
  deletions<-pwalign::deletion(psa)
  for(read_index in seq_along(deletions) ){ # deletions for the reads
    read_del<- unlist(deletions[read_index,])
    for(del_index in seq_along(read_del)){ # for each start location for this read
      psa_mat[read_index, start(read_del[del_index,]):end(read_del[del_index,]) ] <- "Del" # record there's a del for read number i 
    }
  }
  return(psa_mat)
}




