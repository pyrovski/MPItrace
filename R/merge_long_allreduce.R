#!/usr/bin/env Rscript
library(parallel)
library(data.table)

files = Sys.glob('*long_collectives.dat')

f = function(filename){
  a = new.env()
  load(filename, envir=a)
  as.list(a)$s
}

a = mclapply(files, f)
dates = sapply(a, function(x) x$date)
write.table(dates, file='dates.long_collectives.dat', quote=F, row.names=F)
