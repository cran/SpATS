make_design_matrix <-
function(geno, names) {
  Nrow = length(geno)
  Ncol = length(names)
  col = match(geno, names)
  frame = data.frame(i = c(1:Nrow), j = col, v = rep(1,Nrow))
  frame = subset(frame, is.na(col) == FALSE)
  L = as.list(frame)
  X = spam::spam(L, nrow = Nrow, ncol = Ncol)
  return(X)
}
