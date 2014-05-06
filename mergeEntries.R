#!/usr/bin/env Rscript

require('data.table')

readEntry = function(filename){
  e = new.env()
  source(filename, local=e)
  with(e, {
    vars = ls()
    if(length(vars))
      eval(
        parse(
          text=
          paste('data.frame(',
                paste(vars,
                      collapse=','),
                ', stringsAsFactors=F)\n',
                sep='')))
    else NULL
  })
}

if(!exists('f_args')){
  f_args = commandArgs(trailingOnly=T)
  allArgs = commandArgs(trailingOnly=F)
  if(length(f_args) < 2){
    cat(paste("usage:", allArgs[1], "<input file list> <output R file>"))
    quit()
  }
}

inList = readLines(f_args[1])

entries = lapply(inList, readEntry)
classes = unique(rbindlist(lapply(entries,
  function(entry)
  as.data.frame(cbind(name = names(entry),
                      class = sapply(entry, class)),
                stringsAsFactors=F))))
setkey(classes)
entries = rbindlist(lapply(entries, function(entry){
  missing = setdiff(classes$name, names(entry))
  for(col in missing)
##    entry[[col]] = as(NA, classes[J(col)]$class)
    entry[[col]] = NA
  entry = as.data.table(entry)
  setcolorder(entry, classes$name)
  entry
}))

entries[, ranksPerNode:=ceiling(ranks/SLURM_NNODES)]
entries[, SLURM_NNODES:=NULL]

entryCols =
  intersect(c('ranksPerNode','ranks','command'),
            names(entries))
confCols =
  intersect(c('OMP_NUM_THREADS'),
            names(entries))

save(file=f_args[2], entries, entryCols, confCols)
