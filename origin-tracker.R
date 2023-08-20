# !/usr/bin/env Rscript
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("fs", quietly = TRUE)) {
  install.packages("fs", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")
}

library(optparse)
library(fs)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

## parameters
MIN_COUNTRIES_MOB <- 5
DEFAULT_AT_YEAR <- 2021

## read in command-line arguments
option_list = list(
  make_option(c('-d', '--dates'), type='character', default=NULL,
              help='.tsv file containing estimated arrival time for countries of choice, with header (country, date)',
              metavar='character'),
  make_option(c('-o', '--outdir'), type='character', default='./output',
              help='directory where output will be stored', metavar='character'),
  make_option(c('-m', '--mob_data'), type='character', default=sprintf('./resources/shortest_paths_Dmn_matrix_%d.renamed.csv', DEFAULT_AT_YEAR),
              help=sprintf('.csv file containing shortest-path trees D_mn matrix (%d air traffic data is used by default)', DEFAULT_AT_YEAR),
              metavar='character'),
  make_option(c('-n', '--names'), action='store_true', type='logical', default=FALSE,
              help='output names of countries with available mobility data and exit', metavar='logical'),
  make_option(c('-s', '--search'), type='character', default=NULL,
              help='look up countries with name matching input search string and exit', metavar='character'),
  make_option(c('-q', '--query_countries'), type='character', default=NULL,
              help='.txt file containing names of countries for which analysis is to be performed (make sure that countries are separated by newlines, and that there is no header)',
              metavar='character')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## read in shortest-path trees effective distance matrices
cat(sprintf('Reading in shortest-path tree D_mn matrix from %s...', opt$mob_data))
Dmn_mat.df <- read.csv(opt$mob_data)
countries_with_mob <- names(Dmn_mat.df)[2:ncol(Dmn_mat.df)]

## output names of countries with mobility data
if (opt$names) {
  print(countries_with_mob)
  print('Exiting now...')
  quit()
}

## output countries with name matching input search string
if (!is.null(opt$search)) {
  matching_countries <- countries_with_mob[grepl(opt$search, countries_with_mob, ignore.case=TRUE)]
  print(matching_countries)
  print('Exiting now...')
  quit()
}

## read in estimated arrival times for countries of choice from user-specified input file
if (!is.null(opt$dates)) {
  arrival_dates.df <- read.csv(opt$dates, sep='\t')
  cat(sprintf('\nReading in estimated arrival times from %s...', opt$dates))
  
  ## check header consistency
  names(arrival_dates.df) <- tolower(names(arrival_dates.df))
  if (!all(names(arrival_dates.df) == c('country', 'date'))) {
    stop('header in file containing estimated arrival times is invalid, please check')
  }
  
  ## check dates consistency
  arrival_dates.df$date <- as.Date(arrival_dates.df$date)
  if (nrow(arrival_dates.df[is.na(arrival_dates.df$date),]) > 0) {
    stop('estimated arrival dates must be in format YYYY-MM-DD, please check')
  }
  
  ## take only countries with available mobility data
  arrival_dates.mob.df <- arrival_dates.df %>%
    filter(country %in% countries_with_mob)
  ## check if there are at least MIN_COUNTRIES_MOB countries renamining
  if (nrow(arrival_dates.mob.df) < MIN_COUNTRIES_MOB) {
    stop(sprintf('there must be at least %d countries with mobility in order for the analysis to proceed, exiting now...', MIN_COUNTRIES_MOB))
  } else if (nrow(arrival_dates.mob.df) < nrow(arrival_dates.df)) {
    print('WARNING: the following countries are omitted due to absence of mobility data')
    print(arrival_dates.df$country[!arrival_dates.df$country %in% arrival_dates.mob.df$country])
  }
} else {
  stop('.tsv file containing estimated arrival time for countries of choice required')
}

## main analysis
cor_test.p.values <- c()
cor_test.estimates <- c()
mean_Dmns.end <- c()
var_Dmns.end <- c()
cat('\n\nRunning correlation analysis with all countries with mobility data as potential OL, this should take only a couple of seconds...\n')
for (ol in countries_with_mob) {
  ## extract relevant Dmn
  tmp_Dmn_pivoted.df <- Dmn_mat.df[Dmn_mat.df$X == ol,] %>%
    pivot_longer(!X, names_to='country', values_to='Dmn') %>%
    select(-X)
  
  ## construct dataframe
  tmp_Dmn_arrival.compare.df <- arrival_dates.mob.df %>%
    left_join(tmp_Dmn_pivoted.df, by='country') %>%
    mutate(date.diff = as.numeric(
      date - min(arrival_dates.mob.df$date)))
  
  ## correlation test (pearson)
  tmp_cor_test.out <- cor.test(tmp_Dmn_arrival.compare.df[tmp_Dmn_arrival.compare.df$country != ol,]$date.diff,
                               tmp_Dmn_arrival.compare.df[tmp_Dmn_arrival.compare.df$country != ol,]$Dmn,
                               method='pearson')
  
  ## store
  cor_test.p.values <- c(cor_test.p.values, as.numeric(tmp_cor_test.out$p.value))
  cor_test.estimates <- c(cor_test.estimates, as.numeric(tmp_cor_test.out$estimate))
  mean_Dmns.end <- c(mean_Dmns.end, mean(tmp_Dmn_arrival.compare.df$Dmn))
  var_Dmns.end <- c(var_Dmns.end, var(tmp_Dmn_arrival.compare.df$Dmn))
}
## combine as a dataframe
combined_out.df <- data.frame(
  ol_country=countries_with_mob,
  cor_test_p=cor_test.p.values,
  cor_test_est=cor_test.estimates,
  mean_Dmn_end=mean_Dmns.end,
  var_Dmn_end=var_Dmns.end,
  mean_var_end_prod=mean_Dmns.end*var_Dmns.end
)
## sort by cor_test_est
combined_out.df <- combined_out.df %>%
  arrange(-cor_test_est)
combined_out.df$rank <- seq(1:nrow(combined_out.df))

## create output folder if it doesn't already exist
if (dir.exists(opt$outdir)) {
  dir_delete(opt$outdir)
}
dir.create(opt$outdir, showWarnings=TRUE)

## read in names of specific countries to be investigated
query_countries <- c()
if (!is.null(opt$query_countries)) {
  query_countries.df <- read.csv(opt$query_countries, header=FALSE)
  
  ## take only countries with available mobility data
  query_countries.df <- query_countries.df %>%
    filter(V1 %in% countries_with_mob)
  
  ## if no countries remaining, abort
  if (nrow(query_countries.df) > 0) {
    query_countries <- query_countries.df$V1
  } else {
    cat('names of countries specified to be investigated are invalid, please check')
    break
  }
  
  ## create a copy of the dataframe for only query countries
  ## sort by cor_test_est
  query_combined_out.df <- combined_out.df %>%
    filter(ol_country %in% query_countries) %>%
    arrange(-cor_test_est)
  query_combined_out.df$rank <- seq(1:nrow(query_combined_out.df))
  
  ## export output of correlation analysis for query countries to output folder
  cat(sprintf('\nWriting output of correlation analysis for specified query countries to %s/ol-correlation-output.query.csv...', opt$outdir))
  write.csv(query_combined_out.df, file=sprintf('%s/ol-correlation-output.query.csv', opt$outdir), row.names=FALSE, quote=FALSE)
 
  ## plot (query)
  plot.n <- length(query_countries)
  suppressWarnings(g.query <- ggplot() +
                     geom_point(dat=query_combined_out.df,
                                aes(x=rank, y=cor_test_est, size=mean_var_end_prod),
                                shape=21, color='transparent', fill='#c90007', alpha=0.75) +
                     geom_point(dat=query_combined_out.df,
                                aes(x=rank, y=cor_test_est, size=mean_var_end_prod),
                                shape=21, color='#2e2e2e', alpha=0.5) +
                     scale_y_continuous(limits=c(min(query_combined_out.df$cor_test_est)-0.01,
                                                 max(query_combined_out.df$cor_test_est)+0.01)) +
                     scale_x_continuous(expand=c(0, 0), limits=c(-1, plot.n+3)) +
                     coord_cartesian(clip = "off") +
                     labs(x='Rank', y="Correlation coefficient") +
                     theme(
                       legend.position = "none",
                       panel.background = element_rect(fill='transparent'),
                       plot.background = element_rect(fill='transparent', color='NA'),
                       panel.border = element_rect(color='black', size=0.6, fill=NA),
                       panel.grid.major = element_line(color=NA, size=0.),
                       panel.grid.minor = element_line(color=NA, size=0.),
                       # axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.title.y = element_text(size=18, vjust=1),
                       axis.title.x = element_text(size=18)
                     ))
  ## export plot to output folder
  cat(sprintf('\nExporting visualisation of output of correlation analysis for specified query countries to %s/ol-correlation.query.pdf...', opt$outdir))
  ggsave(sprintf('%s/ol-correlation.query.pdf', opt$outdir), device='pdf', width=10, height=6)
}

## export output of correlation analysis to output folder
cat(sprintf('\nWriting output of correlation analysis for specified query countries to %s/ol-correlation-output.global.csv...', opt$outdir))
write.csv(combined_out.df, file=sprintf('%s/ol-correlation-output.global.csv', opt$outdir), row.names=FALSE, quote=FALSE)

## create plot for visualisation
## plot (global)
plot.n <- nrow(combined_out.df)
suppressWarnings(g.global <- ggplot() +
                   geom_point(dat=combined_out.df %>% filter(rank < plot.n),
                              aes(x=rank, y=cor_test_est, size=mean_var_end_prod),
                              shape=21, color='transparent', fill='lightgrey', alpha=0.75) +
                   geom_point(dat=combined_out.df %>% filter(rank < plot.n),
                              aes(x=rank, y=cor_test_est, size=mean_var_end_prod),
                              shape=21, color='#2e2e2e', alpha=0.5) +
                   geom_point(dat=combined_out.df[combined_out.df$ol_country %in% query_countries,],
                              aes(x=rank, y=cor_test_est, size=mean_var_end_prod),
                              shape=21, color='#2e2e2e', fill='#c90007', alpha=0.85) +
                   scale_y_continuous(limits=c(min((combined_out.df %>% filter(rank < plot.n))$cor_test_est)-0.01,
                                               max((combined_out.df %>% filter(rank < plot.n))$cor_test_est)+0.01)) +
                   scale_x_continuous(expand=c(0, 0), limits=c(-3, plot.n+3)) +
                   coord_cartesian(clip = "off") +
                   labs(x='Rank', y="Correlation coefficient") +
                   theme(
                     legend.position = "none",
                     panel.background = element_rect(fill='transparent'),
                     plot.background = element_rect(fill='transparent', color='NA'),
                     panel.border = element_rect(color='black', size=0.6, fill=NA),
                     panel.grid.major = element_line(color=NA, size=0.),
                     panel.grid.minor = element_line(color=NA, size=0.),
                     # axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.y = element_text(size=18, vjust=1),
                     axis.title.x = element_text(size=18)
                   ))
## export plot to output folder
cat(sprintf('\nExporting visualisation of output of correlation analysis for specified query countries to %s/ol-correlation.global.pdf...', opt$outdir))
ggsave(sprintf('%s/ol-correlation.global.pdf', opt$outdir), device='pdf', width=10, height=6)

## exit program
cat('\n\nAll done, exiting now.\n')
