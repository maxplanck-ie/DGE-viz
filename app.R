library(shiny)

# .libPaths('/rstudio/rauer/Rlib_DGE-viz_3.6.0')

# sink('/data/processing4/stdout')
# print(date())
# print(.libPaths())
# print(sessionInfo())
# sink()

lib.dir = '/rstudio/rauer/Rlib_DGE-viz_3.6.0'
library(plotly, lib.loc = lib.dir)
library(data.table, lib.loc = lib.dir)
# library(dplyr, lib.loc = lib.dir)
library(DT, lib.loc = lib.dir)
library(ggplot2, lib.loc = lib.dir)


tab.colnames = NULL
gene_id = NULL

ui <- fluidPage(
  titlePanel("Investigate your RNA-seq DGE data"),
  sidebarLayout(
    sidebarPanel(
      fileInput('inputFile','DGE table'),
      selectInput('inputFormat', label = 'Table format', choices = c('DESeq2','edgeR'), selected = 'DESeq2', multiple = FALSE),
      tabsetPanel(
        tabPanel("Parameter",
                 sliderInput('padj.thrs', 'Threshold for padj', value = 0.05, min = 0, max = .1, step = 0.01),
                 sliderInput('logfc.thrs','Threshold for log2 fold-change', value =c(-1,1),min=-8,max=8,step=1/4),
                 sliderInput('expr.thrs','Thresholds for expression', value=c(-Inf, Inf),min=0, max=15,step=1/10,round=TRUE)
                 # fileInput('annotationFile','Annotation GTF'),
        )
        # tabPanel("Columns",
        #          selectInput('expr_col', label = 'baseMean', choices = tab.colnames),
        #          checkboxInput('logTransform.expr','do log10', value = TRUE),
        #          selectInput('lfc_col', label = 'log2FoldChange', choices = tab.colnames),
        #          checkboxInput('logTransform.foldchange','do log10', value = FALSE),
        #          selectInput('pvalue_col', label = 'pvalue', choices = tab.colnames),
        #          checkboxInput('logTransform.pvalue','do -log10', value = TRUE),
        #          selectInput('padj_col', label = 'padj', choices = tab.colnames),
        #          checkboxInput('logTransform.padj','do -log10', value = TRUE),
        #          actionButton('update_columns',label = "Update")
        # )
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
    print(">>> Loading table")
    tab = fread(input$inputFile$datapath)
    return(tab)
  })
  
  # data_raw = reactive({
  #   tab = fread("/data/iovino/group/snakequest_result/FC_exp9_4_1481_RNAseq_DEseq2_st2_st5_Ac_CG_qnVJmJuH_mRNA-seq/analysis/DESeq2_qnVJmJuH_FC_exp9_4_1481_RNAseq_DEseq2_st2_st5_Ac_CG_st5_HDAC3_sampleSheet/DEseq_basic_DEresults.tsv")
  #   return(tab)
  # })
  
  data_parsed = reactive({
    tab = data_raw()
    print(">>> Parsing table")
    colnames(tab)[1] = 'feature_id'
    if(any(c('external_gene_name', 'symbol') %in% colnames(tab))){
      colnames(tab)[colnames(tab) %in% 'external_gene_name'] = 'symbol'
    } else {
      tab$symbol = tab$feature_id
    }
    cat(">>> Table format:",input$inputFormat,'\n')
    if(!any(c('DESeq2','edgeR') %in% input$inputFormat)) return(NULL)
    
    
    if(input$inputFormat == 'DESeq2'){
      tab.out = data.frame(gene_id = tab$feature_id,
                           baseMean = tab$baseMean, log2FoldChange = tab$log2FoldChange, 
                           pvalue = tab$pvalue, padj = tab$padj,
                           symbol = tab$symbol)
    }
    
    if(input$inputFormat == 'edgeR'){
      tab.out = data.frame(gene_id = tab$feature_id,
                           baseMean = 2**tab$logCPM, log2FoldChange = tab$logFC, 
                           pvalue = tab$PValue, padj = tab$FDR,
                           symbol = tab$symbol)
    }
    
    print('>>> Updating sliders')
    
    updateSliderInput(session, inputId = 'expr.thrs',
                      max = ceiling(max(log2(tab.out$baseMean))),
                      value = c(0,ceiling(max(log2(tab.out$baseMean)))))
    
    lfc.max = ceiling(max(abs(tab.out$log2FoldChange), na.rm = TRUE))
    cat ('>>> lfc slider max:',lfc.max,'\n')
    updateSliderInput(session, inputId = 'logfc.thrs',
                      min = -lfc.max,
                      max = lfc.max)
    
    updateSelectInput(session, inputId = 'genes', choices = setNames(tab.out$gene_id, nm = tab.out$symbol))
    print('>>> Parsing done')
    
    return(tab.out)
  })
  
  data_sliders = reactive({
    tab <- data_parsed()
    print(">>> Parsing selection")
    if(is.null(input$genes)){
      print(">>> Sliders selection") 
      tab0 <- tab %>% mutate(selected = padj < input$padj.thrs & 
                       (log2FoldChange < input$logfc.thrs[1] | log2FoldChange > input$logfc.thrs[2]) &
                       (input$expr.thrs[1] <= log2(baseMean) & log2(baseMean) <= input$expr.thrs[2]))
    } else { 
      print(">>> feature_id selection") 
      tab0 <- tab %>% mutate(selected = (gene_id %in% input$genes) | (symbol %in% input$genes))
    }
    
    print('>>> Content <selected>')
    print(table(tab0$selected, useNA = 'always'))
    
    tab0 <- tab0 %>% mutate(selected = ifelse(is.na(selected), FALSE, selected))
    return(tab0)
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
    tab = data_sliders() 
    if (is.null(d)){
      print(">>> outtab: Slider select")
      tab0 = tab %>% filter(selected) %>% select(-selected)
    } else {
      print(">>> outtab: Key select")
      tab0 = tab %>% filter(gene_id %in% d$key)
    }
    
    formatRound(
      DT::datatable(tab0,
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
      which(sapply(data_parsed(), is.numeric)), digits = 4)
  }) #, server = FALSE)
  
  
  output$preview <- renderTable({
    y = data_raw()[1:5,] %>% as.data.frame
    y = y[,setdiff(colnames(y),'id')]
    y.num = sapply(y, is.numeric)
    y[y.num] = sapply(y[y.num], round, 4)
    return(y)
  })
}

shinyApp(ui, server)
