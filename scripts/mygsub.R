#define mygsub function to replace multiple patterns
mygsub = function(pattern, replacement,x){
  if(length(pattern)!=length(replacement)){
    stop("Pattern and replacement do not have the same length")
  }
  result = x
  for( i in 1:length(pattern)){
    result= gsub(pattern[i],replacement[i],result)
  }
  result
}