library(shiny) # should be provided by shiny server by default
library(DT, lib.loc = './Rlib')
library(ggplot2, lib.loc = './Rlib')
library(data.table, lib.loc = './Rlib')
library(janitor, lib.loc = './Rlib')

source('helpers/helpers.R')
source('helpers/interactiveplots.R')

options(shiny.maxRequestSize=15*1024^2)

tab.colnames = NULL

ui <- fluidPage(
   titlePanel("Investigate your RNA-seq DGE data"),
   sidebarLayout(
     sidebarPanel(
       fileInput('inputFile','DGE table'),
       selectInput('inputFormat', label = 'Table format', choices = c('DEseq2','edgeR','others'), multiple = FALSE),
       tabsetPanel(
         tabPanel("Parameter",
           sliderInput('padj.thrs', 'Threshold for padj', value = 0.05, min = 0, max = 1, step = 1/100),
           sliderInput('logfc.thrs','Threshold for log2 fold-change', value =c(-1,1),min=-8,max=8,step=1/4),
           sliderInput('expr.thrs','Thresholds for expression', value=c(-Inf, Inf),min=0,max=log2(2**13),step=1/10,round=TRUE)
           # fileInput('annotationFile','Annotation GTF'),
         ),
         tabPanel("Columns",
                  selectInput('expr_col', label = 'baseMean', choices = tab.colnames),
                  checkboxInput('logTransform.expr','do log10', value = TRUE),
                  selectInput('lfc_col', label = 'log2FoldChange', choices = tab.colnames),
                  checkboxInput('logTransform.foldchange','do log10', value = FALSE),
                  selectInput('pvalue_col', label = 'pvalue', choices = tab.colnames),
                  checkboxInput('logTransform.pvalue','do -log10', value = TRUE),
                  selectInput('padj_col', label = 'padj', choices = tab.colnames),
                  checkboxInput('logTransform.padj','do -log10', value = TRUE),
                  actionButton('update_columns',label = "Update")
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
       tabPanel("Table - preview input", tableOutput('preview')),
       tabPanel("Sessioninfo", verbatimTextOutput("sessionInfo"))
     )))
)

tablecols.required = c('log2FoldChange','baseMean','pvalue','padj')

server <- function(input, output, session) {
  
  output$sessionInfo <- renderText({
    paste(capture.output(sessionInfo()),collapse = "\n")
  })
  
  data_raw = eventReactive(input$inputFile, { 
    tab = loadData(input$inputFile$datapath)
    tab = clean_names(tab)
    updateSelectInput(session, inputId = 'expr_col', choices = colnames(tab))
    updateSelectInput(session,inputId = 'lfc_col', choices = colnames(tab))
    updateSelectInput(session,inputId = 'pvalue_col', choices = colnames(tab))
    updateSelectInput(session,inputId = 'padj_col', choices = colnames(tab))
    tab$id = paste0('id',nrow(tab))
    tab = tab[order(tab$padj,decreasing = FALSE)[1:10e3], ]
    
    return(tab)
  })
  
  updateColumns <- eventReactive(c(input$inputFormat, input$update_columns),{
    if(input$inputFormat == 'DEseq2')
      return(format.deseq2())
    
    if(input$inputFormat == 'edgeR')
      return(format.edger())
    
    if(input$update_columns){
      cols_idx = sapply(c(input$lfc_col,input$expr_col,input$pvalue_col,input$padj_col),
                        function(x, ref) which(x == ref), colnames(data_raw()))
      print("others - selected columns")
      print(colnames(data_raw())[cols_idx])
      return(cols_idx)
    }
  })

  rows_selected = reactive({
    NULL
  })
  
  data = reactive({
    y = data_raw()
    y = y[,c(1,updateColumns(),grep('id',colnames(y)))]
    colnames(y) <- c('gene_id', tablecols.required, 'id')
    y.copy = y
    
    if(input$logTransform.foldchange)
      y[,tablecols.required[1]] = log10(y[,tablecols.required[1]])
    if(input$logTransform.expr)
      y[,tablecols.required[2]] = log10(y[,tablecols.required[2]])
    if(input$logTransform.pvalue)
      y[,tablecols.required[3]] = -log10(y[,tablecols.required[3]])
    if(input$logTransform.padj)
      y[,tablecols.required[4]] = -log10(y[,tablecols.required[4]])
    
    if(!is.null(input$ma_brush)){
      print("Using ma_brush to select")
      update.vals = ma_brush(input$ma_brush)
      y$selected = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                        update.vals$expr[1] < y$baseMean & y$baseMean < update.vals$expr[2],
                      'selected','not');
    } else if(!is.null(input$volcano_brush)) {
      print("Using volcano_brush to select")
      update.vals = volcano_brush(input$volcano_brush)
      y$selected = ifelse(update.vals$lfc[1] < y$log2FoldChange & y$log2FoldChange < update.vals$lfc[2] &
                            update.vals$pval[1] < y$pvalue & y$pvalue < update.vals$pval[2],
                          'selected','not')

    } else {
      print("Using sliders to select")
      y$selected = ifelse(y.copy$padj < input$padj.thrs &
                            (y$log2FoldChange < input$logfc.thrs[1] | y$log2FoldChange > input$logfc.thrs[2]) &
                            (input$expr.thrs[1] <= y$baseMean & y$baseMean <= input$expr.thrs[2] ),
                          'selected','not');
    }
    y$selected = factor(as.factor(y$selected), levels = c('selected','not',NA));
    y
  })
  
  output$outtab <- renderDT({ 
    x = subset(data(), selected == 'selected')
    x = x[,!(colnames(x) %in% c('selected', 'id'))]
    x.num = sapply(x, class) == class(numeric())
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
  
  output$preview <- renderTable({
    y = data_raw()[1:5,setdiff(colnames(data_raw()), 'id')]
    y.num = sapply(y, is.numeric)
    y[y.num] = apply(y[y.num],2, round, 2)
    head(as.data.frame(y))
  })
  
  output$MAplot <- renderPlot({
    print("Update MA plot")
    dat = data()
    ggplot(dat, aes(baseMean, log2FoldChange)) + geom_point(aes(color = selected)) +
      geom_hline(yintercept = c(input$logfc.thrs), col = 'darkgrey', lty = 2) + 
      geom_vline(xintercept = c(input$expr.thrs), col = 'darkgrey', lty = 2) +
      theme_light()
  })
  
  output$volcanoPlot <- renderPlot({
    print("Update Volcano plot")
    dat = data()
    ggplot(dat, aes(log2FoldChange, pvalue)) + geom_point(aes(color = selected)) +
      geom_vline(xintercept = c(input$logfc.thrs), col = 'darkgrey', lty= 2) + 
      geom_hline(yintercept = min(dat$pvalue[dat$padj < input$padj.thrs], na.rm = TRUE), col = 'darkgrey') + 
      theme_light()
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
