library(readr)
library(reshape2)
library(ggplot2)
library(plotly)

tab0 = read_tsv('./test-data/DEseq_basic_DEresults.fixed.tsv')

fig0 = plot_ly(tab0, x= ~log2(baseMean), y= ~log2FoldChange, 
               color = ~padj < 0.05,
               text = ~paste0(gene_id, ' (',external_gene_name,') <br> -log10 padj:',round(-log10(padj),4)))
fig0

fig1 = plot_ly(tab0, x= ~log2FoldChange, y= ~-log10(pvalue), 
               color = ~padj < 0.05,
               text = ~paste0(gene_id, ' (',external_gene_name,') <br> -log10 padj:',round(-log10(padj),4)))
fig1



p0 = ggplot(tab0, aes(log2(baseMean), log2FoldChange)) + geom_point(aes(color = padj < 0.05)) + 
  theme_classic()
fig0.a = ggplotly(p0)
toWebGL(p = fig0.a)


# useful collection

# Source: https://plotly-r.com/client-side-linking.html
# - highlight on click
highlight(on = "plotly_hover", off = "plotly_doubleclick")

# add select by text field
highlight(
  # time_series, 
  on = "plotly_click", 
  selectize = TRUE, # add text field
  dynamic = TRUE, 
  persistent = TRUE
)

