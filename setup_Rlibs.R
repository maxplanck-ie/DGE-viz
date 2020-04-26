rversion = gsub('R version (.*) (.*)','\\1',sessionInfo()$R.version$version.string)
lib.dir = paste0('Rlib_',rversion)

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

dir.create(lib.dir, showWarnings = TRUE)

b = !(cran.packages %in% rownames(installed.packages(lib.loc = lib.dir)))

if(all(!b))
  stop('All cran.packages installed already. Exit')

cat(paste0(">> Installing packages to \'",lib.dir,"\'\n"), cran.packages[b])
install.packages(cran.packages[b], lib = lib.dir, Ncpus = 4)

# for(pkg0 in cran.packages[b])
#   cat(paste0('install.packages(\"',pkg0,'\", lib = \'./Rlib/\', Ncpus = 4)'),'\n')
