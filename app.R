library(shiny)
library(DT)
library(ggplot2)

source('helpers/helpers.R')
source('helpers/interactiveplots.R')

ui <- fluidPage(
   titlePanel("Investigate your RNA-seq DGE data"),
   sidebarLayout(
     sidebarPanel(
       fileInput('inputFile','DGE table'),
       tabsetPanel(
         tabPanel("Parameter",
           sliderInput('padj.thrs', 'Threshold for padj', value =0.05, min = 0, max = 1, step = 1/100),
           sliderInput('logfc.thrs','Threshold for log2 fold-change', value =c(-1,1),min=-8,max=8,step=1/4),
           sliderInput('expr.thrs','Thresholds for expression', value=c(-Inf, Inf),min=0,max=log2(2**13),step=1/10,round=TRUE)
           # fileInput('annotationFile','Annotation GTF'),
         ),
         tabPanel("Columns",
                  h1("Warning:"),
                  p("Change options only if your table is not canonical DESeq2 format"),
                  numericInput('expr_col', label = 'baseMean', value = 2, min = 1),
                  numericInput('lfc_col', label = 'log2FoldChange', value = 3, min = 1),
                  numericInput('pvalue_col', label = 'pvalue', value = 6, min = 1),
                  numericInput('padj_col', label = 'padj', value = 7, min = 1),
                  actionButton('updateidx','Update')
                )
       ),
       verbatimTextOutput('errorOut', placeholder = TRUE)
   ),
   mainPanel(
     tabsetPanel(
       tabPanel('Plots',
                plotOutput("MAplot", brush = "ma_brush"),
                plotOutput("volcanoPlot", brush = "volcano_brush")),
       tabPanel('Table - selected data', DT::dataTableOutput('outtab')),
       tabPanel("Table - preview input", tableOutput('preview'))
     )))
)

tablecols.required = c('log2FoldChange','baseMean','pvalue','padj')

server <- function(input, output, session) {
  
  data_raw = eventReactive(input$inputFile, { 
    loadData(input$inputFile$datapath)
  })
  
  data = reactive({
      y = data_raw()
      
      idx.cols = c(input$lfc_col,input$expr_col,input$pvalue_col,input$padj_col)
      colnames(y)[c(input$lfc_col,input$expr_col,input$pvalue_col,input$padj_col)] <- tablecols.required
      
      if(!is.null(input$ma_brush)){
        print("Using ma_brush to select")
        update.vals = ma_brush(input$ma_brush)
        y$sign = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                          2**update.vals$expr[1] < y$baseMean & y$baseMean < 2**update.vals$expr[2],
                        'significant','not');
      } else if(!is.null(input$volcano_brush)) {
        print("Using volcano_brush to select")
        update.vals = volcano_brush(input$volcano_brush)
        y$sign = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                          update.vals$pval[1] < -log10(y$pvalue) & -log10(y$pvalue) < update.vals$pval[2],
                        'significant','not');
      } else {
        print("Using sliders to select")
        y$sign = ifelse(y$padj < input$padj.thrs & 
                          ( y$log2FoldChange < input$logfc.thrs[1] | y$log2FoldChange > input$logfc.thrs[2]) &
                          (2**input$expr.thrs[1] <= y$baseMean & y$baseMean <= 2**input$expr.thrs[2] ),
                        'significant','not');
      }
      y$sign = factor(as.factor(y$sign), levels = c('significant','not',NA));
      y
  })
  
  output$preview <- renderTable({
    y = data_raw()[1:5,]
    y.num = apply(y, 2, is.numeric)
    y[,y.num] = apply(y[,y.num],2, round, 2)
    head(y)
  })
  
  output$MAplot <- renderPlot({
    print("Update MA plot")
    dat = data()
    ggplot(dat, aes(log2(baseMean), log2FoldChange)) + geom_point(aes(color = sign)) +
      geom_hline(yintercept = c(input$logfc.thrs), col = 'darkgrey', lty = 2) + 
      geom_vline(xintercept = c(input$expr.thrs), col = 'darkgrey', lty = 2) +
      theme_light()
  })
  
  output$volcanoPlot <- renderPlot({
    print("Update Volcano plot")
    dat = data()
    ggplot(dat, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(color = sign)) +
      geom_vline(xintercept = c(input$logfc.thrs), col = 'darkgrey', lty= 2) + 
      geom_hline(yintercept = min(-log10(dat$pvalue[dat$padj < input$padj.thrs]), na.rm = TRUE), col = 'darkgrey') + 
      theme_light()
  })
  
  output$outtab <- renderDT({ 
    x = subset(data(), sign == 'significant');
    
    x = x[,!(colnames(x) %in% 'sign')]
    
    x.num = sapply(dd1, class) == class(numeric())
    # x[,x.num] = apply(x[,x.num],2,round,4)
    
    
    formatRound(
      DT::datatable(x,
                    rownames = FALSE,
                    extensions = 'Buttons', options = list(
                      dom = 'Bfrtip',
                      buttons = list('copy', 
                                     list(                        
                                       extend = 'collection',
                                       buttons = c('csv', 'excel', 'pdf'),
                                       text = 'Download'
                                     ))
                    )),
                which(x.num), digits = 2)
    }, server=FALSE)
  output$errorOut <- renderText({
    b1 = tablecols.required %in% colnames(data_raw())
    if(!all(b1))
      return(paste("Columns not found:\n", paste(tablecols.required[!b1], sep = '\n')))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
