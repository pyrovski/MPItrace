#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')
source('./read.R')

require('data.table')

## group experiments by entryCols from entries

load('mergedEntries.Rsave')

confSpace = unique(entries[,entryCols,with=F])
setkey(confSpace)
setkeyv(entries, entryCols)

countedConfSpace = entries[confSpace, list(count=nrow(.SD))]

sel = which(!complete.cases(confSpace))
if(length(sel)){
  cat('Removing', length(sel), 'cases:\n')
  print(countedConfSpace[sel])
}
confSpace = na.omit(confSpace)

confSpace$key = rowApply(confSpace, toKey)
setkeyv(confSpace, entryCols)

countedConfSpace = entries[confSpace, list(count=nrow(.SD))]

g = function(entry){
  filename = file.path(entry$path, 'merged.Rsave')
  ## Does merged data exist?  If not, merge it.
  if(-1 == file.access(filename)){
    print(filename)
    run(path=entry$path, saveResult=T)
  } ##else cat('Merged data for', entry$date, 'already exists\n')
  
  ## Read merged data.
  e = new.env()
  load(filename, envir=e)
  result = as.list(e)
  result$date = entry$date
  return(result)
}

f = function(conf){
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

### Merge runtimes tables
  runtimes = rbindlist(lapply(result, function(r){
    r$runtimes$date = r$date
    r$runtimes
  }))

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
  setkey(runtimes, date, rank, uid)
  ranks = (1:conf$ranks)-1
  hashMaps = nnapply(ranks, function(r){
    result = as.data.table(nnapply(dates, function(d){
      unique(runtimes[J(d, r)]$hash)
    }))
  })

  fixDates = tail(dates, -1)
  masterDate = head(dates, 1)
  if(length(fixDates) &&
     any(sapply(hashMaps, function(h)
                any(unlist(rowApply(h, function(row)
                                    length(unique(unlist(row)))>1)))))){
    ##!@todo test on runs with different hashes
    cat('Matching hashes between', length(dates), 'runs\n')
    setkey(runtimes, date, rank, hash)
    lapply(ranks, function(r){
      hashes = hashMaps[[as.character(r)]]
      masterHashes = hashes[[masterDate]]
      lapply(fixDates, function(d){
        key = data.table(date = d, rank = r, hash = hashes[[d]])
        masterKey = data.table(date = masterDate, rank = r, hash = masterHashes)
        setkey(key)
        setkey(masterKey)
        runtimes[key]$hash = runtimes[masterKey]$hash
      })
    })
    cat('Done matching hashes\n')
  }

### Match requests between runs
  
  runtimes
}

##!@todo this may run into memory limitations. If so, just run the
##!configurations one at a time.

go = function(){
  merged <<- mcrowApply(confSpace, f)
  names(merged) <<- confSpace$key
  save(merged, file='mergedData.Rsave')
}

if(!interactive())
  go()
