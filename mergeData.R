#!/usr/bin/env Rscript

##!@todo group by contents of command file
options(mc.cores=16)

source('~/local/bin/pbutils.R')
source('./read.R')
source('./util.R')

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
    result$messageEdges[[col]] = entry[[col]]
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

  confName = gsub('[/.]', '_', conf$key)
  write.table(unique(runtimes[,confCols,with=F]),
              file=paste(confName,'confSpace.csv',sep='.'), quote=F,
              sep=',', row.names=F)

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
  ##x$messageEdges = x$messageEdges[[1]]
  x$messageEdges = rbindlist(x$messageEdges)
  by = c('s_uid','d_uid',confCols)
  x$messageEdges = x$messageEdges[,lapply(.SD, mean),by=by]

  x$compEdges = rbindlist(x$compEdges)
  by = c('s_uid',confCols)
  x$compEdges = x$compEdges[,lapply(.SD, mean),by=by]
  setkey(x$messageEdges)
  x$edges = merge(x$messageEdges, x$compEdges, all=T)

  ## set power to zero for message edges
  x$edges[is.na(power), power:=0]
  
  x$vertices = x$reduced[!is.na(name),list(name=head(name,1)),by=vertex]
  
  ## Get an initial schedule, starting with minimum time per task.
  x$schedule = x$edges[,.SD[which.min(weight)],by=list(s_uid)]

  schedule = getSchedule(x$schedule)
  g = schedule$g
  x$schedule = schedule$edges
  rm(schedule)
  
  setkey(x$schedule, start, weight)
  x$schedule$e_uid = 1:nrow(x$schedule)

  ## critical path
  setkey(x$schedule, src)

  ## store e_uid for edges on the longest path.
  ## edges in topological order: x$schedule[J(x$vertices)]
  x$longestPath = longest.path(x$schedule[J(x$vertices)], x$vertices, g)
  rm(g)

  ## copy e_uid to edges
  setkey(x$schedule, s_uid, d_uid)
  setkey(x$edges, s_uid, d_uid)
  x$edges[x$schedule, e_uid:=e_uid]

  ## Insert slack edges. These edges will have power, but not minimum
  ## time. To insert the new edges, we need new vertices and new edge
  ## uids. We use the negative of the original edge uid for each.
  x$edges = slackEdges(x$edges, x$longestPath)
  x$edges[is.na(weight), weight:=0]

  ## recompute schedule with slack edges
  x$schedule = x$edges[,.SD[which.min(weight)],by=list(e_uid)]
  ## get topological order again and add start times
  schedule = getSchedule(x$schedule)
  x$schedule = schedule$edges
  x$vertices = schedule$vertices
  rm(schedule)

  ##!@todo get a src and dest rank for each edge to facilitate slack
  ##!attribution. This should happen in edges, not schedule.
  
  ## setkey(x$reduced, uid)
  ## x$schedule$s_rank = x$reduced[J(x$schedule[, s_uid]), rank, mult='first']
  ## x$schedule$d_rank = x$reduced[J(x$schedule[, d_uid]), rank, mult='first']

  
  confName = gsub('[/.]', '_', x$key)
  ## write edge uids
  write.table(unique(x$edges[,list(e_uid)])[order(e_uid)],
              file=paste(confName,'task_IDs.csv',sep='.'),
              row.names=F, quote=F, sep=',')
  
  write.table(x$vertices,
              file=paste(confName, '.vertices.csv', sep=''),
              row.names=F, quote=F, sep=',')

  firstCols = c('e_uid', confCols)
  setcolorder(x$edges, c(firstCols, setdiff(names(x$edges), firstCols)))
  
  write.table(x$edges[,c(firstCols, 'src', 'dest', 'weight', 'power'),with=F],
              file=paste(confName, '.edges.csv', sep=''),
              row.names=F, quote=F, sep=',')

  ##!@todo I would like to keep all the edges in one table for the LP,
  ##!but this may be impossible.
  
  ## write.table(x$edges[,c(firstCols, 'src', 'dest', 'weight', 'power'),with=F],
  ##             file=paste(confName, '.slack_edges.csv', sep=''),
  ##             row.names=F, quote=F, sep=',')
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

  confSpace <<- unique(entries[,confCols,with=F])
  countedEntryspace <<- entries[entrySpace, list(count=nrow(.SD)),
                                by=entryCols]

  sel = which(!complete.cases(entrySpace))
  if(length(sel)){
    cat('Removing', length(sel), 'cases:\n')
    print(countedEntryspace[sel])
  }
  entrySpace <<- na.omit(entrySpace)

  entrySpace$key <<- rowApply(entrySpace, toKey)
  setkeyv(entrySpace, entryCols)

  countedEntryspace <<- entries[entrySpace, list(count=nrow(.SD)), by=entryCols]
  measurementCols <<- c('duration','pkg_w','pp0_w','dram_w')

  cat('Merging configurations\n')
  merged <<- mcrowApply(entrySpace, mergeConfs)
  names(merged) <<- entrySpace$key
  lapply(names(merged), function(name) merged[[name]]$key<<-name)
  cat('Done merging configurations\n')
  cat('Reducing configurations\n')
  reduced <<- mclapply(merged, reduceConfs)
  cat('Done reducing configurations\n')
  save(measurementCols, reduced, merged, entrySpace, countedEntryspace,
       entryCols, entries, confSpace, confCols,
       file='mergedData.Rsave')
}

if(!interactive())
  go()
