#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')
require('data.table')
require('parallel')

if(!exists('f_args') && !interactive()){
  f_args = commandArgs(trailingOnly=T)
  allArgs = commandArgs(trailingOnly=F)
  if(length(f_args) < 2){
    cat(paste("usage:", allArgs[1], "<input file list> <output R file>\n"))
    quit()
  }
}

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

mergeEntries = function(inList = readLines(f_args[1]), outFile = f_args[2]){
  entries <<- mclapply(inList, readEntry)
  errors <<- sapply(entries, is.null)
  #errorFiles <<- inList[errors]
  entries <<- entries[!errors]
#  lengths <<- sapply(entries, length)
#  errors <<- lengths < max(lengths)
#  entries <<- entries[!errors]
  classes <<- unique(rbindlist(lapply(entries,
    function(entry)
    as.data.frame(cbind(name = names(entry),
                        class = sapply(entry, class)),
                  stringsAsFactors=F))))
  setkey(classes)
  entries <<- rbindlist(mclapply(entries, function(entry){
    missing = setdiff(classes$name, names(entry))
    for(col in missing)
      ##    entry[[col]] = as(NA, classes[J(col)]$class)
      entry[[col]] = as(NA, classes[J(col), class])
    entry = as.data.table(entry)
    setcolorder(entry, classes$name)
    entry
  }))
  
  entries[, ranksPerNode:=ceiling(ranks/SLURM_NNODES)]
  entries[, SLURM_NNODES:=NULL]
  
  entryCols <<-
    intersect(c('ranksPerNode','ranks','command'),
              names(entries))
  confCols <<-
    intersect(
      c('OMP_NUM_THREADS', ## number of OpenMP threads
        'cpuFreq'          ## static CPU frequency in kHz
        ),
      names(entries))

  entrySpace <<- unique(entries[,entryCols,with=F])

  setkey(entrySpace)
  setkeyv(entries, entryCols)

  countedEntryspace <<- entries[entrySpace, list(count=nrow(.SD)),
                                by=entryCols]

  sel = which(!complete.cases(entrySpace))
  if(length(sel)){
    cat('Removing', length(sel), 'cases:\n')
    print(countedEntryspace[sel])
  }
  entrySpace <<- na.omit(entrySpace)

  entrySpace$key <<- unlist(rowApply(entrySpace, toKey))
  setkeyv(entrySpace, entryCols)

  countedEntryspace <<-
    entries[entrySpace, list(count=nrow(.SD)), by=entryCols]

  save(file=outFile, entries, entryCols, confCols, entrySpace,
       countedEntryspace)
}

if(!interactive())
  mergeEntries()
##mergeEntries(readLines('entries'), 'mergedEntries.Rsave')
