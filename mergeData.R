#!/usr/bin/env Rscript

##!@todo group by contents of command file
options(mc.cores=16)

source('~/local/bin/pbutils.R')
source('./read.R')
source('./util.R')

require('data.table')

getEntryData = function(entry){
  filename = file.path(entry$path, 'merged.Rsave')
  ## Does merged data exist?  If not, merge it.
  if(-1 == file.access(filename)){
    print(filename)
    ##! if on lc system, categorize.sh can launch moab jobs for this
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

mergeConfs = function(conf, entries){
  print(conf)
  result = rowApply(entries[conf], getEntryData)
  
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
  cat('Reducing between configs\n')
  if(T){
    nonMeasurementCols = setdiff(names(x$runtimes), measurementCols)
    n_time = system.time(
      x$reduced <-
      x$runtimes[,c(measurementCols, by), with=F][,
                                            lapply(.SD[,measurementCols,
                                                       with=F], mean),
                                            by=by]
      )
    cat('reduce time: ', n_time[3], '\n')
    x$runtimes = x$runtimes[,nonMeasurementCols,with=F]
    setkeyv(x$reduced, by)
    setkeyv(x$runtimes, by)
    x$reduced = x$runtimes[x$reduced,,mult='first']
    n_time = system.time(
      x$reduced <-
      reduceNoEffect(x$reduced, measurementCols, nonMeasurementCols, by))
    cat('reduceNoEffect time: ', n_time[3], '\n')
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
 
  cat('Message edges\n')
  ## all message edges should be identical between runs
  ##x$messageEdges = x$messageEdges[[1]]
  x$messageEdges = rbindlist(x$messageEdges)
  by = c('s_uid','d_uid',confCols)
  x$messageEdges = x$messageEdges[,lapply(.SD, mean),by=by]
  x$messageEdges[, type:='message']
  setkey(x$reduced, uid)
  ## message edges can convert to slack edges, so we need the destination rank
  x$messageEdges[, rank:=x$reduced[J(x$messageEdges[,d_uid,by=d_uid]), rank]]
  x$messageEdges =
    reduceNoEffect(x$messageEdges, c('weight'),
                   setdiff(names(x$messageEdges),
                           c('weight')), c('s_uid','d_uid',confCols))
  
  cat('Computation edges\n')
  x$compEdges = rbindlist(x$compEdges)
  by = c('s_uid',confCols)
  x$compEdges = x$compEdges[,lapply(.SD, mean),by=by]
  x$compEdges[, type:='comp']
  x$compEdges[, rank:=x$reduced[J(x$compEdges[,d_uid,by=d_uid]), rank]]
  x$compEdges =
    reduceNoEffect(x$compEdges, c('weight','power'),
                   setdiff(names(x$compEdges),
                           c('weight','power')), c('s_uid','d_uid',confCols))

  setkey(x$messageEdges)
  setkey(x$compEdges, NULL)
  x$edges = merge(x$messageEdges, x$compEdges, all=T)

  ## set power to zero for message edges
  x$edges[is.na(power), power:=0]

  x$vertices = x$reduced[!is.na(name),list(name=head(name,1)),by=vertex]
  
  ## Get an initial schedule, starting with minimum time per task.
  x$schedule = x$edges[,.SD[which.min(weight)],by=list(s_uid)]

  cat('Schedule and critical path\n')
  schedule = getSchedule(x$schedule)
  ##g = schedule$g
  x$schedule = schedule$edges
  x$critPath = schedule$critPath
  rm(schedule)
  
  ## copy e_uid to edges
  setkey(x$schedule, s_uid, d_uid, type)
  setkey(x$edges, s_uid, d_uid, type)
  x$edges[x$schedule, e_uid:=e_uid]

  cat('Pareto frontiers\n')
  ## get pareto frontiers
  x$edges = pareto(x$edges)

  ## Insert slack edges. These edges will have power, but not minimum
  ## time. To insert the new edges, we need new vertices and new edge
  ## uids. We use the negative of the original edge uid for each.
  cat('Slack edges\n')
  x$edges = slackEdges(x$edges, x$critPath)
  x$edges[is.na(weight), weight:=0]

  ## Recompute schedule with slack edges. This will not change the
  ## critical path because no slack edges are on it.
  x$schedule = x$edges[,.SD[which.min(weight)],by=list(e_uid)]
  ## get topological order again and add start times
  schedule = getSchedule(x$schedule, doCritPath=F)
  x$schedule = schedule$edges
  x$vertices = schedule$vertices
  rm(schedule)
  return(x)
}

writeSlice = function(x){
  x$edges = timeslice(x$schedule, x$edges)[[1]]
  
  firstCols = c('e_uid', confCols)
  setcolorder(x$edges, c(firstCols, setdiff(names(x$edges), firstCols)))
  confName = gsub('[/.]', '_', x$key)

###!@todo for now, only output one or two configs per edge. This
###!allows a convex representation of the time-power relationship over
###!the configuration space. We already output only the configs on the
###!pareto frontier; in the future, let the optimizer handle a
###!piecewise-linear time-power function.
  
  ## write confSpace
  write.table(unique(x$edges[,confCols,with=F], by=confCols),
              file=paste(confName,'confSpace.csv',sep='.'), quote=F,
              sep=',', row.names=F)

  x$vertices = x$edges[,list(vertex=union(src,dest))]
  
  write.table(x$vertices,
              file=paste(confName, '.vertices.csv', sep=''),
              row.names=F, quote=F, sep=',')

###!@todo for dependent timeslices, only output one configuration for
###!the first computation edge on each rank. This configuration should
###!come from the previous timeslice solution.
  
  write.table(x$edges[,list(src=head(src, 1), dest=head(dest, 1),
                            edge_rank=head(rank,1)),
                      by=e_uid],
              file=paste(confName, '.edges.csv', sep=''),
              row.names=F, quote=F, sep=',')
  write.table(x$edges[,c(firstCols, 'weight', 'power'),with=F],
              file=paste(confName, '.edge_weights.csv', sep=''),
              row.names=F, quote=F, sep=',')

  ## ##!I only need a single rank column. Slack edges always go on the
  ## ##!destination rank, and we don't care about messages. Only the
  ## ##!computation and slack edges should have rank information.

  ##!@todo this might choose the wrong edge for a rank if multiple
  ##!edges start at the same time? I guess we should break ties by
  ##!choosing computation or slack edges, but not message
  ##!edges. Alternatively, the timeslice time can just depend on all
  ##!tasks.
  write.table(
    x$schedule[e_uid %in% x$edges[,e_uid]][type %in% c('comp', 'slack'),
                                           .SD[which.max(start)],
                                           by=rank][,list(rank, last_edge=e_uid)],
    file=paste(confName, '.last_edges.csv', sep=''),
    row.names=F, quote=F, sep=',')
}

##!@todo this may run into memory limitations. If so, just run the
##!configurations one at a time.

go = function(){
  ## group experiments by entryCols from entries

  load('mergedEntries.Rsave', envir=.GlobalEnv)

  ##!
  entries = entries[ranks==8]

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

  countedEntryspace <<-
    entries[entrySpace, list(count=nrow(.SD)), by=entryCols]
  measurementCols <<- c('duration','pkg_w','pp0_w','dram_w')

  f = function(entry){
    cat('Merging configurations\n')
    merged <- mergeConfs(entry, entries)
    cat('Done merging configurations\n')
    cat('Reducing configurations\n')
    reduced <- reduceConfs(merged)
    cat('Done reducing configurations\n')
    cat('Writing timeslices\n')
    writeSlice(reduced)
    cat('Done writing timeslices\n')
    return(list(merged=merged, reduced=reduced))
  }
  result = mcrowApply(entrySpace, f)
  
  
  save(measurementCols, result, entrySpace, countedEntryspace,
       entryCols, entries, confSpace, confCols,
       file='mergedData.Rsave')
}

if(!interactive())
  go()
