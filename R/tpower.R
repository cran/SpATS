tpower <-
function(x, t, p) {
# Function for truncated p-th power function
   return((x - t) ^ p * (x > t))
}
