#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')
source('./read.R')

require('data.table')

## group experiments by confCols from entries

load('mergedEntries.Rsave')

confSpace = unique(entries[,confCols,with=F])
setkey(confSpace)
setkeyv(entries, confCols)

countedConfSpace = entries[confSpace, list(count=nrow(.SD)), by=confCols]

sel = which(!complete.cases(confSpace))
if(length(sel)){
  cat('Removing', length(sel), 'cases:\n')
  print(countedConfSpace[sel])
}
confSpace = na.omit(confSpace)

confSpace$key = rowApply(confSpace, toKey)

g = function(entry){
  filename = file.path(entry$path, 'merged.Rsave')
  if(-1 == file.access(filename)){
    print(filename)
    run(path=entry$path, saveResult=T)
  } ##else cat('Merged data for', entry$date, 'already exists\n')
  
  ## Read merged data.
  e = new.env()
  load(filename, envir=e)
  as.list(e)
}

f = function(conf){
  ## Does merged data exist?  If not, merge it.
  print(conf)
  rowApply(entries[conf], g)
}
merged = mcrowApply(confSpace, f)
names(merged) = confSpace$key
