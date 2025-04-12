randomly_flip_edges <- function(adj_matrix, flip_prop = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Ensure no self-loops
  diag(adj_matrix) <- 0
  
  # Find existing edges in the upper triangle
  upper_indices <- which(upper.tri(adj_matrix) & adj_matrix == 1, arr.ind = TRUE)
  
  # Determine number of edges to flip
  num_to_flip <- ceiling(flip_prop * nrow(upper_indices))
  if (num_to_flip == 0) return(adj_matrix)  # nothing to flip
  
  # Randomly sample indices to flip
  flip_indices <- upper_indices[sample(1:nrow(upper_indices), num_to_flip), , drop = FALSE]
  
  # Flip selected edges to -1 symmetrically
  for (i in seq_len(nrow(flip_indices))) {
    row <- flip_indices[i, 1]
    col <- flip_indices[i, 2]
    adj_matrix[row, col] <- -1
    adj_matrix[col, row] <- -1
  }
  
  return(adj_matrix)
}
