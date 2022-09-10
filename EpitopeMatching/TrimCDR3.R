# Trim a CDR3 sequence. Only trim the first character if it's a "C". Only trim
# the last character if the first character is a "C" AND the last character is 
# an "F" or a "W". This has to be done separately for chain 1 (alpha) and chain 
# 2 (beta).
TrimCDR3 <- function(cdr3) {
  front <- lapply(cdr3, str_sub, 0, 1)
  back <- lapply(cdr3, str_sub, -1)
  
  trim1 <- which(front == "C")
  trim2 <- intersect(which(back == "F" | back == "W"), trim1)
  
  cdr3.trimmed <- cdr3
  cdr3.trimmed[trim1] <- str_sub(cdr3.trimmed[trim1], 2)
  cdr3.trimmed[trim2] <- str_sub(cdr3.trimmed[trim2], 0, -2)
  return(cdr3.trimmed)
}