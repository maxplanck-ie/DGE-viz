parse_yaml <- function(config_file){
  t0 = read.table(config_file, sep = ':', header = FALSE, stringsAsFactors = FALSE)
  v0 = setNames(gsub(' +','',t0[,2]), nm = t0[,1])
  v0
}
