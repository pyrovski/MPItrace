#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')
require('data.table')
require('parallel')

## if(!exists('f_args') && !interactive()){
##   f_args = commandArgs(trailingOnly=T)
##   allArgs = commandArgs(trailingOnly=F)
##   if(length(f_args) < 1){
##     cat(paste("usage:", allArgs[1], "<input list file>\n"))
##     quit(status=1)
##   }
## }
## print(f_args)

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

mergeEntries = function(inList = readLines(f_args[1]), outFile = 'mergedEntries.Rsave'){
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

  ## add gmpi_replay, gmpi_replay_file, powerLimit, (powerBalancing?) to classes
  cols =
    c('ranksPerNode','ranks','command',
###!@todo this assumes gmpi_replay_file is only set if gmpi_replay is set
      'gmpi_replay','gmpi_replay_file',
      'powerLimit'
      #, 'powerBalancing'
      )
  colClasses = list(ranksPerNode='numeric', ranks='numeric', command='character',
    gmpi_replay='numeric', gmpi_replay_file='character', powerLimit='numeric')
    #powerBalancing='numeric')
  for(col in cols)
    if(!col %in% classes$name)
      classes = rbind(classes, data.table(name=col, class=colClasses[[col]]))
  
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
    intersect(cols,
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

  ## countedEntrySpace <<- entries[entrySpace, list(count=nrow(.SD)),
  ##                               by=entryCols]

  ## sel = which(!complete.cases(entrySpace))
  ## if(length(sel)){
  ##   cat('Removing', length(sel), 'cases:\n')
  ##   print(countedEntrySpace[sel])
  ## }
  ## entrySpace <<- na.omit(entrySpace)

  entrySpace$key <<- unlist(rowApply(entrySpace, toKey))
  setkeyv(entrySpace, entryCols)

  countedEntrySpace <<-
    entries[entrySpace, list(count=nrow(.SD)), by=entryCols]

  save(file=outFile, entries, entryCols, confCols, entrySpace,
       countedEntrySpace)
}

#if(!interactive())
  mergeEntries(readLines('entries'))
