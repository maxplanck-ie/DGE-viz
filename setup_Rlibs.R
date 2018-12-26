packages = c('DT',
             'data.table',
             'dplyr',
             'readr',
             'ggplot2',
             'rlang',
             'crosstalk',
             'htmlwidgets',
             'tidyr',
             'tidyselect',
             'tidyverse')

dir.create("./", showWarnings = TRUE)

b = !(packages %in% rownames(installed.packages(lib.loc = './Rlib')))

if(all(!b))
  stop('All packages installed already. Exit')


for(p in packages[b]){
  cat(">> Installing",p,' to ./Rlib\n')
  install.packages(p, lib = './Rlib/')
}
