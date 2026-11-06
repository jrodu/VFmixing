# The estimation functions print progress bars and the selected model.
# Run an expression without that output and return its value.
quietly = function(expr){
  res = NULL
  utils::capture.output(res <- expr)
  res
}
