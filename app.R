library(shiny)
library(data.table)
library(dplyr)
library(plotly)
library(DT)

tab.colnames = NULL
gene_id = NULL

ui <- fluidPage(
  titlePanel("Investigate your RNA-seq DGE data"),
  sidebarLayout(
    sidebarPanel(
      fileInput('inputFile','DGE table'),
      selectInput('inputFormat', label = 'Table format', choices = c('DEseq2','edgeR','others'), multiple = FALSE),
      tabsetPanel(
        tabPanel("Parameter",
                 sliderInput('padj.thrs', 'Threshold for padj', value = 0.05, min = 0, max = .1, step = 0.01),
                 sliderInput('logfc.thrs','Threshold for log2 fold-change', value =c(-1,1),min=-8,max=8,step=1/4),
                 sliderInput('expr.thrs','Thresholds for expression', value=c(-Inf, Inf),min=0,max=log2(2**15),step=1/10,round=TRUE)
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
      selectizeInput(
        inputId = "genes", 
        label = NULL,
        # placeholder is enabled when 1st choice is an empty string 
        choices = c("Choose a gene" = "", gene_id),
        multiple = TRUE
      ),
      tabsetPanel(
        tabPanel('Plots',
                 plotlyOutput(outputId = "ma"),
                 plotlyOutput(outputId = "volcano")),
        tabPanel('Table - selected data', DT::dataTableOutput('outtab')),
        tabPanel("Table - preview input", tableOutput('preview'))
        # ,
        # tabPanel("Sessioninfo", verbatimTextOutput("sessionInfo"))
      )))
)

server <- function(input, output, session, ...) {
  
  data_raw = eventReactive(input$inputFile, { 
    tab = fread(input$inputFile$datapath)
    if(colnames(tab)[1] == 'V1')
      colnames(tab)[1] = 'gene_id'
    
    
    # tab = clean_names(tab)
    updateSelectInput(session, inputId = 'expr_col', choices = colnames(tab))
    updateSelectInput(session, inputId = 'lfc_col', choices = colnames(tab))
    updateSelectInput(session, inputId = 'pvalue_col', choices = colnames(tab))
    updateSelectInput(session, inputId = 'padj_col', choices = colnames(tab))
    
    updateSliderInput(session, inputId = 'expr.thrs',
                      max = ceiling(max(log2(tab$baseMean))),
                      value = c(0,ceiling(max(log2(tab$baseMean)))))
    updateSelectInput(session, inputId = 'genes', choices = setNames(tab$gene_id, nm = tab$symbol))
    
    return(tab)
  })
  
  data_sliders = reactive({
    if(is.null(input$genes))
      return(data_raw() %>% mutate(selected = padj < input$padj.thrs & 
                              (log2FoldChange < input$logfc.thrs[1] | log2FoldChange > input$logfc.thrs[2]) &
                              (input$expr.thrs[1] <= log2(baseMean) & log2(baseMean) <= input$expr.thrs[2])))
    
    return(data_raw() %>% mutate(selected = gene_id %in% input$genes | symbol %in% input$genes))
  })
  
  output$ma <- renderPlotly({
    
    p1 <- ggplot(data = data_sliders(),  aes(log2(baseMean), log2FoldChange)) + 
      geom_point(aes(color = selected,
                     text = paste0(symbol,' (',gene_id,')'),
                     key = gene_id), show.legend = FALSE) + 
      scale_color_manual(values = c('TRUE' = 'blue','FALSE'='grey')) + 
      geom_hline(yintercept = input$logfc.thrs, color = 'darkgrey', lty = 2) +
      geom_vline(xintercept = input$expr.thrs, color = 'darkgrey', lty = 2) + 
      theme_light()
    
    height <- session$clientData$output_p_height
    width <- session$clientData$output_p_width
    
    toWebGL(ggplotly(p1, height = height, width = width, tooltip = c('text','x','y')))
  })
  
  output$volcano <- renderPlotly({
    
    p2 <- ggplot(data = data_sliders(), aes(log2FoldChange, -log10(pvalue))) +
      geom_point(aes(color = selected, 
                     text = paste0(symbol,' (',gene_id,')'),
                     key = gene_id), show.legend = FALSE) + 
      scale_color_manual(values = c('TRUE' = 'blue','FALSE'='grey')) + 
      geom_hline(yintercept = input$padj.thrs, color = 'darkgrey', lty = 2) + 
      geom_vline(xintercept = input$logfc.thrs, color = 'darkgrey', lty = 2) +
      theme_light()
    
    height <- session$clientData$output_p_height
    width <- session$clientData$output_p_width
    
    toWebGL(ggplotly(p2, height = height, width = width, tooltip = c('text','x','y')))
  })  
  
  # output$p <- renderPlotly({
  #   # req(input$genes)
  #   # if (identical(input$genes, "")) return(NULL)
  #   
  #   p1 <- ggplot(data = tab0, aes(log2(baseMean), log2FoldChange)) + 
  #     geom_point(aes(color = external_gene_name %in% input$genes, 
  #                    text = paste0(external_gene_name,' (',gene_id,')'),
  #                    key = gene_id)) + 
  #     scale_color_manual(values = c('TRUE' = 'blue','FALSE'='grey'))
  #   p2 <- ggplot(data = tab0, aes(log2FoldChange, -log10(pvalue))) +
  #     geom_point(aes(color = external_gene_name %in% input$genes, 
  #                    text = paste0(external_gene_name,' (',gene_id,')'),
  #                    key = gene_id)) + 
  #     scale_color_manual(values = c('TRUE' = 'blue','FALSE'='grey'))
  #   height <- session$clientData$output_p_height
  #   width <- session$clientData$output_p_width
  #   
  #   toWebGL(subplot(which_layout = 'merge', 
  #                   ggplotly(p1, height = height, width = width, tooltip = c('text','x','y')),
  #                   ggplotly(p2, height = height, width = width, tooltip = c('text','x','y')),
  #           nrows = 2))
  # })
  
  output$outtab <- renderDataTable({
    d <- event_data("plotly_selected")
    if (is.null(d))
      tab0 = data_sliders() %>% select(-selected)
    else 
      tab0 = data_raw() %>% filter(gene_id %in% d$key)
    formatRound(
      DT::datatable(tab0 ,
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
      which(sapply(data_raw(), is.numeric)), digits = 4)
  }, server = FALSE)
  
  
  output$preview <- renderTable({
    print(data_raw())
    y = data_raw()[1:5,] %>% as.data.frame
    y = y[,setdiff(colnames(y),'id')]
    y.num = sapply(y, is.numeric)
    print(y.num)
    y[y.num] = sapply(y[y.num], round, 4)
    return(y)
  })
}

shinyApp(ui, server)