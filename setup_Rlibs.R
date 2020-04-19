cran.packages = c('DT',
             'data.table',
             'dplyr',
             'readr',
             'ggplot2',
             'rlang',
             'crosstalk',
             'htmlwidgets',
             'tidyr',
             'tidyselect',
             'janitor',
             'snakecase')

dir.create("./", showWarnings = TRUE)

b = !(cran.packages %in% rownames(installed.packages(lib.loc = './Rlib')))

if(all(!b))
  stop('All cran.packages installed already. Exit')


for(p in cran.packages[b]){
  cat(">> Installing",p,' to ./Rlib\n')
  install.packages(p, lib = './Rlib/')
}
