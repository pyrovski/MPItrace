#!/usr/bin/env Rscript

readRuntime = function(filename){
  a = data.table(read.table(filename,h=T,stringsAsFactors=F))
  a[size == -1,]$size = NA
  a[dest == -1,]$dest = NA
  a[src == -1,]$src = NA
  a[comm == 0,]$comm = NA
  ##!@todo convert request list from string to R list
}
