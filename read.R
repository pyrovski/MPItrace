#!/usr/bin/env Rscript

library('data.table')
source('~/local/bin/pbutils.R')

readRuntime = function(filename){
  firstLine = strsplit(readLines(filename, n=1),'\t')[[1]]
  colClasses = c(
    start='numeric',
    duration='numeric',
    name='character',
    size='integer',
    dest='integer',
    src='integer',
    tag='integer',
    comm='character',
    hash='character',
    pkg_w='numeric',
    pp0_w='numeric',
    dram_w='numeric',
    reqs='character')
  a = data.table(read.table(filename,h=T,stringsAsFactors=F,colClasses=colClasses))
  a[size == -1,]$size = NA
  a[dest == -1,]$dest = NA
  a[src == -1,]$src = NA
  a[is.na(size)][tag == -1]$tag = NA
  a$reqs = strsplit(a$reqs,',')
  return(a)
}

readAll = function(path='.'){
  files = sort(list.files(path,'runtime.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  result = lapply(files, readRuntime)
  names(result) = ranks
  result = napply(result, function(x, name){x$rank = as.numeric(name);x})
  result = rbindlist(result)
  return(result)
}

deps = function(x){
  
}
