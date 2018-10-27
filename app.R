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
       sliderInput('padj.thrs', 'Threshold for padj', value =0.05, min = 0, max = 1, step = 1/100),
       sliderInput('logfc.thrs','Threshold for log2 fold-change', value =c(-1,1),min=-8,max=8,step=1/4),
       sliderInput('expr.thrs','Thresholds for expression', value=c(-Inf, Inf),min=0,max=log2(2**13),step=1/10,round=TRUE),
       # fileInput('annotationFile','Annotation GTF'),
       verbatimTextOutput('errorOut')
     ),
     mainPanel(
       tabsetPanel(
         tabPanel('Plots',
                  plotOutput("MAplot", brush = "ma_brush"),
                  plotOutput("volcanoPlot", brush = "volcano_brush")),
         tabPanel('Table', DT::dataTableOutput('outtab'))
     )))
)

tablecols.required = c('log2FoldChange','baseMean','pvalue','padj')

server <- function(input, output, session) {
  data_raw = eventReactive( input$inputFile, { 
    loadData(input$inputFile$datapath)
  })
  
  data = reactive({
      y = data_raw()
      if(!is.null(input$ma_brush)){
        print("Using ma_brush to select")
        update.vals = ma_brush(input$ma_brush)
        y$sign = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                          2**update.vals$expr[1] < y$baseMean & y$baseMean < 2**update.vals$expr[2],
                        'significant','not');
        y$sign = relevel(as.factor(y$sign), ref = 'significant');
      } else if(!is.null(input$volcano_brush)) {
        print("Using volcano_brush to select")
        update.vals = volcano_brush(input$volcano_brush)
        y$sign = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                          update.vals$pval[1] < -log10(y$pvalue) & -log10(y$pvalue) < update.vals$pval[2],
                        'significant','not');
        y$sign = relevel(as.factor(y$sign), ref = 'significant');
      } else {
        print("Using sliders to select")
        y$sign = ifelse(y$padj < input$padj.thrs & 
                          ( y$log2FoldChange < input$logfc.thrs[1] | y$log2FoldChange > input$logfc.thrs[2]) &
                          (2**input$expr.thrs[1] <= y$baseMean & y$baseMean <= 2**input$expr.thrs[2] ),
                        'significant','not');
        y$sign = relevel(as.factor(y$sign), ref = 'significant');
      }
      y
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
    DT::datatable(round(x[order(x$padj, decreasing = FALSE),], 4), 
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                  ))
  }, server=FALSE)
  
  output$errorOut <- renderText({
    b1 = tablecols.required %in% colnames(data_raw())
    if(!all(b1))
      return(paste("Columns not found:\n", paste(tablecols.required[!b1], sep = '\n')))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
