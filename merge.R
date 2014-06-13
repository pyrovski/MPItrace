#!/usr/bin/env Rscript
#msub -l nodes=1

options(mc.cores=16)
args = commandArgs(trailingOnly=T)
if(length(args) < 1){
  stop('expected path argument')
}
paths = args

source('read.R')

lapply(paths, function(path){
  tryCatch(run(path=path, saveResult=T, noReturn=T),
           error = function(e) e)
})
