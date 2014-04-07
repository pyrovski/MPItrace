#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')

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

f = function(row){
  ## Does merged data exist?
    
  ## If not, merge it.

  ## Read merged data.
}
merged = rowApply(confSpace, f)
names(merged) = confSpace$key
