xy_range_str <- function(e) {
  if(is.null(e)) return("NULL\n")
  paste0("xmin=", round(e$xmin, 2), " xmax=", round(e$xmax, 2), 
         " ymin=", round(e$ymin, 2), " ymax=", round(e$ymax, 2))
}

ma_brush <- function(e){
  print('interactive ma_brush')
  if(is.null(e)) return(list())
  return(list(expr = c(round(e$xmin, 2), round(e$xmax, 2)), 
              lfc = c(round(e$ymin, 2), round(e$ymax, 2))))
}

volcano_brush <- function(e){
  if(is.null(e)) return(list())
  return(list(lfc = c(round(e$xmin, 2), round(e$xmax, 2)),
              pval = c(round(e$ymin, 2), round(e$ymax, 2))))
}
