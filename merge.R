#!/usr/bin/env Rscript
#msub -l nodes=1

args = commandArgs(trailingOnly=T)
if(length(args) != 1){
  stop('expected path argument')
}
path = args[[1]]

source('read.R')

setwd(path)
run(saveResult=T)
