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

# tablecols.required = c('log2FoldChange','baseMean','pvalue','padj')

format.deseq2 <- function(){
  # DESeq2 table - baseMean, log2FC, pvalue,padj
  return(c(3,2,6,7))
}

format.edger <- function(){
  # edger table - baseMean, log2FC, pvalue,padj
  return(c(3,2,5,6))
}