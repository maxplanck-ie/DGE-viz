source('../helpers/helpers.R')
config0 = parse_yaml('../config.yaml')
# BiocManager::install("pasilla", lib = '/rstudio/rauer/Rlib_DGE-viz_3.6.0')
# BiocManager::install("DESeq2", lib = '/rstudio/rauer/Rlib_DGE-viz_3.6.0')
# BiocManager::install("edgeR", lib = '/rstudio/rauer/Rlib_DGE-viz_3.6.0')

lib.dir0 = config0['lib.path']
if(is.na(lib.dir0))
  lib.dir0 = NULL
library(DESeq2, lib.loc = lib.dir0)
library(pasilla, lib.loc = lib.dir0)
library(dplyr, lib.loc = lib.dir0)
library(tibble, lib.loc = lib.dir0)
library(readr, lib.loc = lib.dir0)
library(edgeR, lib.loc = lib.dir0)

pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds = DESeq(dds)

tab0 = results(dds) %>% as.data.frame() %>% rownames_to_column()
write_tsv(tab0,'./Pasilla_testdata.DESeq2.tsv')

tab1 = lfcShrink(dds,coef="condition_untreated_vs_treated", type="normal") %>% as.data.frame() %>% rownames_to_column()
write_tsv(tab0,'./Pasilla_testdata.DESeq2_lfcShrink.tsv')

# edgeR
y <- DGEList(counts=cts,group=as.factor(coldata$type))
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~as.factor(coldata$type))
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tab2 = topTags(qlf, n = Inf) %>% as.data.frame %>% rownames_to_column()
write_tsv(tab2,'./Pasilla_testdata.edgeR.tsv')
