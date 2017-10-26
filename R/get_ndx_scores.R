get_ndx_scores <-
function(weights, geno, names) {
  scored = as.vector(weights != 0)
  geno2 = subset(geno, scored == TRUE)
  col = match(geno2,names)
  u = as.vector(na.omit(unique(col)))
  return (sort(u))
}
