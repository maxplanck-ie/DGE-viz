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

dir.create("./Rlib", showWarnings = TRUE)

b = !(cran.packages %in% rownames(installed.packages(lib.loc = './Rlib')))

if(all(!b))
  stop('All cran.packages installed already. Exit')

cat(">> Installing packages to ./Rlib\n", cran.packages[b])
install.packages(cran.packages[b], lib = './Rlib/', libs_only = TRUE)
