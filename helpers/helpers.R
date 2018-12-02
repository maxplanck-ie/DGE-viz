parseData <- function(data.file){
  require(data.table, lib.loc = './Rlib')
  return(as.data.frame(fread(data.file)))
}

loadSimulatedData <- function(){
  source('./helpers/simulator.R')
  deseq2.tab = produceData()
}

loadData <- function(tab.file = NULL){
  if(is.null(tab.file)){
    return(loadSimulatedData())
  } else {
    return(parseData(tab.file))
  }
}
