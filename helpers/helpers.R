
parseData <- function(data.file){
  require(tidyverse)
  return(read_tsv(data.file))
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

xy_range_str <- function(e) {
  if(is.null(e)) return("NULL\n")
  paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
         " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
}