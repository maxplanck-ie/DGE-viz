parseData <- function(data.file){
  require(data.table, lib.loc = './Rlib')
  return(as.data.frame(fread(data.file)))
}

loadData <- function(tab.file){
  return(parseData(tab.file))
}
