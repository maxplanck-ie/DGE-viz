
produceData <- function(){
  require(DESeq2)
  dd1 = makeExampleDESeqDataSet(betaSD = 2)
  keep = rowSums(counts(dd1)) > 10
  dd1 = dd1[keep,]
  
  dd1 = DESeq(dd1)
  res = results(dd1)
  return(as.data.frame(res))  
}
   