#!/usr/bin/env Rscript

##!@todo group by contents of command file
options(mc.cores=16)

source('~/local/bin/pbutils.R')
source('./read.R')

require('data.table')

g = function(entry){
  filename = file.path(entry$path, 'merged.Rsave')
  ## Does merged data exist?  If not, merge it.
  if(-1 == file.access(filename)){
    print(filename)
    ##!@todo if on lc system, launch moab jobs for this
    run(path=entry$path, saveResult=T)
  } ##else cat('Merged data for', entry$date, 'already exists\n')
  
  ## Read merged data.
  e = new.env()
  tryCatch(load(filename, envir=e),
           error=function(e) cat('fin:', filename, entry$date, str(e), '\n'),
           finally=NULL)
  result = as.list(e)
  result$date = entry$date
  for(col in confCols){
    result$runtimes[[col]] = entry[[col]]
    result$compEdges[[col]] = entry[[col]]
  }
  return(result)
}

mergeConfs = function(conf){
  print(conf)
  result = rowApply(entries[conf], g)
  
  dates = sapply(result, '[[', 'date')

### Merge between runs of the same config. This entails matching
### hashes, requests, and communicators between runs. After the
### matching, the runs should have identical per-rank event
### ordering. UIDs will not match between runs unless we force them to
### match.

  ## Sanity test
  if(length(unique(sapply(result, function(r) nrow(r$runtimes)))) > 1){
    errMsg = 'Runs differ in number of events!'
    cat(errMsg, '\n')
    stop(errMsg)
  }

  f = function(r, name){
    r[[name]]$date = r$date
    r[[name]]
  }
  
### Merge runtimes tables
  runtimes = rbindlist(lapply(result, f, name='runtimes'))

  assignments = rbindlist(lapply(result, f, name='assignments'))

  messageEdges = lapply(result, '[[', 'messageEdges')
  names(messageEdges) = lapply(result, '[[', 'date')
  compEdges = lapply(result, '[[', 'compEdges')
  names(compEdges) = lapply(result, '[[', 'date')

  ##!@todo merge runtimes by config, then recreate graphs from merged runtimes
  
  rm(result)
 
### Comms should already be unified; I mapped MPI_COMM_WORLD and
### MPI_COMM_NULL to -1 and -2, respectively. The other comms should
### align, but we need to test this.
  setkey(runtimes, date, rank)
  commSeqs = lapply(dates, function(d){
    runtimes[J(d)][comm != '(nil)']$comm
  })
  if(length(unique(commSeqs)) != 1){
    errMsg = "FIXME match comms between runs"
    cat(errMsg, '\n')
    stop(errMsg)
  }
    
### Construct hash maps; one table per rank with one column per run
  ## setkey(runtimes, date, rank, uid)
  ## ranks = (1:conf$ranks)-1
  ## hashMaps = nnapply(ranks, function(r){
  ##   result = as.data.table(nnapply(dates, function(d){
  ##     unique(runtimes[J(d, r)]$hash)
  ##   }))
  ## })

  fixDates = tail(dates, -1)
  masterDate = head(dates, 1)
  cat('Matching hashes between', length(dates), 'runs\n')
  hashMaps = runtimes[date == masterDate, list(uid, rank, hash)]
  setkey(hashMaps)
  runtimes$hash = NULL
  setkey(runtimes, uid, rank)
  runtimes = runtimes[hashMaps]
  cat('Done matching hashes\n')

### Compute power consumption
  
### Match requests between runs?

###!@todo match UIDs between runs
  ## uidCheck =
  ##   runtimes[,list(name, rank, uid, hash)][,lapply(.SD,function(col)
  ##                                                  length(unique(col))),
  ##                                          .SDcols=c('name','hash'),
  ##                                          by=list(uid)]
  ## if(any(uidCheck[,2:ncol(uidCheck),with=F] != 1)){
  ##   stop(conf$key, ' failed UID check')
  ## }
  
  list(runtimes=runtimes, assignments=assignments,
       messageEdges=messageEdges, compEdges=compEdges)
}

## combine within confCols combinations. This will combine multiple
## runs with the same configuration.
reduceConfs = function(x){
  x$runtimes$date = NULL
  x$runtimes$start = NULL
  by = c('uid',confCols)
  if(T){
    nonMeasurementCols = setdiff(names(x$runtimes), measurementCols)
    x$reduced =
      x$runtimes[,lapply(.SD[,measurementCols,
                             with=F], mean),
                 by=by]
    x$runtimes = x$runtimes[,nonMeasurementCols,with=F]
    setkeyv(x$reduced, by)
    setkeyv(x$runtimes, by)
    x$reduced = x$runtimes[x$reduced,,mult='first']
  } else {
    ## this is slower by 3x.
    nonMeasurementCols = setdiff(names(x$runtimes), c(by, measurementCols))
    x$reduced =
      x$runtimes[,c(lapply(.SD[,measurementCols,with=F],
                           mean),
                    lapply(.SD[,nonMeasurementCols,with=F],
                           function(col) head(col,1))),
                 by=by]
  }
  x$runtimes = NULL
  ## all message edges should be identical between runs
  x$messageEdges = x$messageEdges[[1]]
  x$compEdges = rbindlist(x$compEdges)
  x$compEdges = x$compEdges[,lapply(.SD, mean),by=by]
  setkey(x$messageEdges)
  x$edges = merge(x$messageEdges, x$compEdges, all=T)
  return(x)
}

##!@todo this may run into memory limitations. If so, just run the
##!configurations one at a time.

go = function(){
  ## group experiments by entryCols from entries

  load('mergedEntries.Rsave', envir=.GlobalEnv)

  entrySpace <<- unique(entries[,entryCols,with=F])
  setkey(entrySpace)
  setkeyv(entries, entryCols)

  countedEntryspace <<- entries[entrySpace, list(count=nrow(.SD)), by=entryCols]

  sel = which(!complete.cases(entrySpace))
  if(length(sel)){
    cat('Removing', length(sel), 'cases:\n')
    print(countedEntryspace[sel])
  }
  entrySpace <<- na.omit(entrySpace)

  entrySpace$key = rowApply(entrySpace, toKey)
  setkeyv(entrySpace, entryCols)

  countedEntryspace <<- entries[entrySpace, list(count=nrow(.SD)), by=entryCols]
  measurementCols <<- c('duration','pkg_w','pp0_w','dram_w')

  cat('Merging configurations\n')
  merged <<- mcrowApply(entrySpace, mergeConfs)
  names(merged) <<- entrySpace$key
  cat('Done merging configurations\n')
  cat('Reducing configurations\n')
  reduced <<- mclapply(merged, reduceConfs)
  cat('Done reducing configurations\n')
  save(measurementCols, reduced, merged, entrySpace, countedEntryspace,
       entryCols, entries, file='mergedData.Rsave')
}

if(!interactive())
  go()
