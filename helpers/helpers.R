parseData <- function(data.file){
  require(data.table, lib.loc = './Rlib')
  
  return(as.data.frame(fread(data.file)))
}

loadData <- function(tab.file){
  return(parseData(tab.file))
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
