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
  result$compEdges[, c('src', 'dest') :=
                   list(as.integer(src), as.integer(dest))]
  ##!@todo set messageEdges to NULL in read.R if no messages
  if(!inherits(result$messageEdges, 'logical'))
    result$messageEdges[, c('src', 'dest') :=
                        list(as.integer(src), as.integer(dest))]
  for(col in confCols){
    result$compEdges[[col]] = entry[[col]]
    ##! message edges should be unaffected. If actual message times
    ##are ever measured, this should change.
    ##if(!any(is.na(result$messageEdges)))
    ##  result$messageEdges[[col]] = entry[[col]]
  }
  return(result)
}

mergeConfs = function(conf, entries){
  print(conf)
  setkeyv(entries, entryCols)
  cat(length(entries[conf, which=T]), ' entries\n')
  result = rowApply(entries[conf], getEntryData)

  cat(conf$key, 'loaded data\n')
  
  dates = sapply(result, '[[', 'date')

### Merge between runs of the same config. This entails matching
### hashes, requests, and communicators between runs. After the
### matching, the runs should have identical per-rank event
### ordering. UIDs will not match between runs unless we force them to
### match.

  ## Sanity test
  if(length(unique(sapply(result, function(r) nrow(r$compEdges)))) > 1){
    errMsg = 'Runs differ in number of edges!'
    cat(errMsg, '\n')
    stop(errMsg)
  }

  f = function(r, name){
    r[[name]]$date = r$date
    r[[name]]
  }
  
  assignments = rbindlist(lapply(result, f, name='assignments'))
  for(i in 1:length(result)) result[[i]]$assignments=NULL
  gc()

  messageEdges = lapply(result, '[[', 'messageEdges')
  names(messageEdges) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$messageEdges=NULL
  gc()

  compEdges = lapply(result, '[[', 'compEdges')
  names(compEdges) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$compEdges=NULL
  gc()

  collectives = lapply(result, '[[', 'collectives')
  names(compEdges) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$collectives=NULL
  gc()
  
  rm(result)
 
### Comms should already be unified; I mapped MPI_COMM_WORLD and
### MPI_COMM_NULL to -1 and -2, respectively. The other comms should
### align, but we need to test this.

###!@todo I disabled these because we don't have enough RAM to load
###!the runtimes table for all runs in a sequence.
  ## setkey(runtimes, date, rank)
  ## commSeqs = lapply(dates, function(d){
  ##   runtimes[J(d)][comm != '(nil)']$comm
  ## })
  ## if(length(unique(commSeqs)) != 1){
  ##   errMsg = "FIXME match comms between runs"
  ##   cat(errMsg, '\n')
  ##   stop(errMsg)
  ## }
    
  ## fixDates = tail(dates, -1)
  ## masterDate = head(dates, 1)
  ## cat('Matching hashes between', length(dates), 'runs\n')
  ## hashMaps = runtimes[date == masterDate, list(uid, rank, hash)]
  ## setkey(hashMaps)
  ## runtimes$hash = NULL
  ## setkey(runtimes, uid, rank)
  ## runtimes = runtimes[hashMaps]
  ## cat('Done matching hashes\n')
 
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

  list(##runtimes=runtimes,
       assignments=assignments,
       messageEdges=messageEdges,
       compEdges=compEdges,
       collectives=collectives)
}

## combine within confCols combinations. This will combine multiple
## runs with the same configuration.
reduceConfs = function(x){
  startTime = Sys.time()

  by = c('uid',confCols)
  cat('Reducing between configs\n')

  x$collectives = unique(x$collectives)
  if(length(x$collectives) > 1){
    stop('Unmatched collectives\n')
  }
  x$collectives = x$collectives[[1]]

  ## all message edges should be identical between runs
  if(any(is.na(x$messageEdges))){
    cat('No message edges in at least one run\n')
    x$messageEdges = NULL
  } else {
    cat('Merging message edges\n')
    measurementCols = c('weight')
### this unique is ok because message weights are identical between
### runs
    x$messageEdges = unique(rbindlist(x$messageEdges))
    cat('Message edges:', nrow(x$messageEdges), '\n')
    by = c('s_uid','d_uid')
    cores = getOption('mc.cores')
    if(!is.null(cores) && cores > 1){
      setkey(x$messageEdges, s_uid)
      s_uids = unique(x$messageEdges[, s_uid])
      s_uid_chunks = chunk(s_uids, length(s_uids)/cores)
      rm(s_uids)
      x$messageEdges = rbindlist(mclapply(s_uid_chunks, function(s_uids)
        x$messageEdges[J(s_uids)][,lapply(.SD, mean),by=by]))
      rm(s_uid_chunks)
    } else
      x$messageEdges = x$messageEdges[,lapply(.SD, mean),by=by]
    x$messageEdges_inv =
      x$messageEdges[, setdiff(names(x$messageEdges), measurementCols), with=F]
    x$messageEdges = x$messageEdges[, c('s_uid', 'd_uid', measurementCols), with=F]
  }
  
  cat('Computation edges:', sum(sapply(x$compEdges, nrow)), '\n')
  measurementCols = c('weight','power')
  nonMeasurementCols =
    setdiff(names(x$compEdges[[1]]), c(confCols,measurementCols))
  ##!@todo this assumes all compEdges have the same structure
  x$compEdges_inv =
    x$compEdges[[1]][, nonMeasurementCols, with=F][, head(.SD, 1), by=s_uid]
  setkey(x$compEdges_inv, s_uid)

  ##!@todo this could be done in read.R

  ##!@todo d_uid is redundant for computation edges, but useful for
  ##merging with message edges
  
  by = c('s_uid','d_uid',confCols)
  x$compEdges =
    rbindlist(mclapply(x$compEdges, function(compEdges)
                       compEdges[, c(by, measurementCols), with=F]
                       ))
  setkey(x$compEdges, s_uid)
  cores = getOption('mc.cores')
  if(!is.null(cores) && cores > 1){
    s_uids = unique(x$compEdges[, s_uid])
    s_uid_chunks = chunk(s_uids, length(s_uids)/cores)
    rm(s_uids)
    x$compEdges = rbindlist(mclapply(s_uid_chunks, function(s_uids)
      x$compEdges[J(s_uids)][,lapply(.SD, mean),by=by]))
    rm(s_uid_chunks)
  } else
    x$compEdges = x$compEdges[,lapply(.SD, mean),by=by]
  setkey(x$compEdges, s_uid)
  x$compEdges =
    reduceNoEffect(x$compEdges, x$compEdges_inv, measurementCols,
                   by, 's_uid')
  x$compEdges_inv[, type:='comp']

  if(!is.null(x$messageEdges)){
    x$messageEdges_inv[, type:='message']
    commonNames = intersect(names(x$messageEdges_inv), names(x$compEdges_inv))
    messageOnlyNames = setdiff(commonNames, names(x$messageEdges_inv))
    compOnlyNames = setdiff(commonNames, names(x$compEdges_inv))
    setkeyv(x$messageEdges_inv, commonNames)
    setkey(x$compEdges_inv, NULL)
    x$edges_inv = merge(x$messageEdges_inv, x$compEdges_inv, all=T)
  } else
    x$edges_inv = x$compEdges_inv
  x$compEdges_inv = NULL
  x$messageEdges_inv = NULL

  ## assign edge uids
  if(!'e_uid' %in% names(x$edges_inv)){
    setkey(x$edges_inv, s_uid, d_uid, type)
    uids = x$edges_inv[, list(e_uid=.GRP), keyby=list(s_uid, d_uid, type)]
    e = x$edges_inv[uids]
    x$edges_inv = e
    rm(uids, e)
    gc()
  }
  ##! replace s_uid/d_uid with s_uid in _inv tables as key
  setkey(x$edges_inv, s_uid, d_uid)
  setkey(x$compEdges, s_uid, d_uid)
  x$compEdges = x$compEdges[x$edges_inv[,list(s_uid, d_uid, e_uid)]]
  x$compEdges[, c('s_uid','d_uid') := list(NULL, NULL)]
  setkey(x$compEdges, e_uid)

  if(!is.null(x$messageEdges)){
    setkey(x$messageEdges, s_uid, d_uid)
    x$messageEdges = x$messageEdges[x$edges_inv[,list(s_uid, d_uid, e_uid)]]  
    x$messageEdges[, c('s_uid','d_uid') := list(NULL, NULL)]
    setkey(x$messageEdges, e_uid)
  }
  setkey(x$edges_inv, e_uid)
  
  cat('Pareto frontiers\n')
  ## get pareto frontiers
  x$compEdges = pareto(x$compEdges)
  gc()

###!@todo split here; save and restore to get around parallel memory
###issue
  
### merge x$compEdges and x$messageEdges; this involves using extra
### space for message edge configs, but makes edge lookup easier
  if(!is.null(x$messageEdges)){
    messageMinConf = minConf(x$compEdges)[, confCols, with=F]
    x$messageEdges[, confCols := messageMinConf, with=F]
    x$messageEdges[, power := 0]
    x$edges = rbind(x$compEdges, x$messageEdges, use.names=T)
  } else
    x$edges = x$compEdges
  x$compEdges = NULL
  x$messageEdges = NULL
  gc()

###! renumber vertices from 1:n
  x$vertices = x$edges_inv[, list(vertex=union(src, dest))][order(vertex)]
  x$vertices$newVertex = 1:nrow(x$vertices)
  setkey(x$vertices, vertex)
  setkey(x$edges_inv, src)
  x$edges_inv = x$edges_inv[x$vertices]
  x$edges_inv = x$edges_inv[, c('src', 'newVertex') := list(newVertex, NULL)]
  setkey(x$edges_inv, dest)
  x$edges_inv = x$edges_inv[x$vertices]
  x$edges_inv = x$edges_inv[, c('dest', 'newVertex') := list(newVertex, NULL)]
  x$vertices[, c('vertex', 'newVertex') := list(newVertex, NULL)]
  
  ## Get an initial schedule, starting with minimum time per task.
  x$schedule = x$edges[,.SD[which.min(weight)],keyby=e_uid]
  ## get src and dest vertices
  x$schedule = x$schedule[x$edges_inv[, list(e_uid, src, dest, rank, type)]]
  
  cat('Schedule and critical path\n')
  schedule = getSchedule(x$schedule, x$vertices)
  ##g = schedule$g
  x$schedule = schedule$edges
  x$critPath = schedule$critPath
  x$vertices = schedule$vertices
  rm(schedule); gc()

  ## Insert slack edges. These edges will have power, but not minimum
  ## time. To insert the new edges, we need new vertices and new edge
  ## uids. We use the negative of the original edge uid for each.
  cat('Slack edges\n')
  x$schedule = x$schedule[x$edges_inv[, list(e_uid, s_uid, d_uid)]]
  x$schedule = slackEdges(x$schedule, x$critPath)
  x$schedule[is.na(weight), weight:=0]

  ## add slack vertices to vertices table
  x$slackVertices =
    x$schedule[is.na(s_uid),
               list(start=head(start, 1),via=-head(src, 1)), by=src]
  setnames(x$slackVertices, c('vertex', 'start', 'via'))
  return(x)
}

writeSlices = function(x, sliceDir='csv'){
  dir.create(sliceDir)
  options(scipen=7)
  firstCols = c('e_uid', confCols)
  confName = gsub('[/.]', '_', x$key)

###!@todo this would use less memory if we computed one slice at a time
  slices = timeslice(x$schedule, rbind(x$vertices, x$slackVertices), x$edges)
  names(slices) = sprintf('%.3f', as.numeric(names(slices)))
  writeSlice = function(slice, sliceTime){
    sliceName = paste(confName, sliceTime, sep='_')
    setcolorder(slice, c(firstCols, setdiff(names(slice), firstCols)))

###!@todo for now, only output one or two configs per edge. This
###!allows a convex representation of the time-power relationship over
###!the configuration space. We already output only the configs on the
###!pareto frontier; in the future, let the optimizer handle a
###!piecewise-linear time-power function.
    
    ## write confSpace
    write.table(unique(slice[,confCols,with=F], by=confCols),
                file=file.path(sliceDir, paste(sliceName,'confSpace.csv',sep='.')),
                quote=F, sep=',', row.names=F)

    write.table(slice[,list(vertex=union(src,dest))],
                file=file.path(sliceDir, paste(sliceName, '.vertices.csv', sep='')),
                row.names=F, quote=F, sep=',')

###!@todo for dependent timeslices, only output one configuration for
###!the first computation edge on each rank. This configuration should
###!come from the previous timeslice solution.

    write.table(slice[,list(
      ##!@todo this is easier to get from x$schedule
      src=head(src, 1), dest=head(dest, 1),
      edge_rank=head(rank,1),
      ## but not these
      minTime=min(weight),
      maxTime=max(weight),
      minPower=min(power),
      maxPower=max(power)),
                      by=e_uid],
                file=file.path(sliceDir, paste(sliceName, '.edges.csv', sep='')),
                row.names=F, quote=F, sep=',')
    write.table(slice[,c(firstCols, 'weight', 'power'),with=F],
                file=file.path(sliceDir, paste(sliceName, '.edge_weights.csv', sep='')),
                row.names=F, quote=F, sep=',')
    write.table(slice[,list(count=nrow(.SD)), by=e_uid][count > 1, list(e_uid)],
                file=file.path(sliceDir, paste(sliceName, '.edge_multiConf.csv', sep='')),
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
      x$schedule[e_uid %in% slice[,e_uid]][type %in% c('comp', 'slack'),
                                           .SD[which.max(start)],
                                           by=rank][,list(rank, last_edge=e_uid)],
      file=file.path(sliceDir, paste(sliceName, '.last_edges.csv', sep='')),
      row.names=F, quote=F, sep=',')
  }
  mclapply(names(slices), function(sliceTime) writeSlice(slices[[sliceTime]], sliceTime))
  NULL
}

##!@todo this may run into memory limitations. If so, just run the
##!configurations one at a time.

go = function(){
  ## group experiments by entryCols from entries

  load('mergedEntries.Rsave', envir=.GlobalEnv)

  measurementCols <<- c('duration','pkg_w','pp0_w','dram_w')
  confSpace <<- unique(entries[,confCols,with=F])
  setkey(confSpace)

  f = function(entry){
    startTime = Sys.time()
    filename = paste('mergedData', gsub('[/.]', '_', entry$key),
      'Rsave', sep='.')    
    if(file.exists(filename)){
      cat(entry$key, 'already merged\n')
      return(NULL)
    }

### check for missing configs
    confs = unique(entries[entry, confCols, with=F])
    setkey(confs)
    missingConfs = confSpace[!confs]
    if(nrow(missingConfs)){
      warning(cat(entry$key, 'missing', nrow(missingConfs), 'configs\n'),
              immediate.=T)
      print(missingConfs)
      save(missingConfs, file=paste('missing', gsub('[/.]', '_', entry$key), sep='.'))
    }
    
    cat(entry$key, 'Merging configurations\n')
    merged <- mergeConfs(entry, entries)
    cat(entry$key, 'Done merging configurations\n')
    cat(entry$key, 'merge time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Reducing configurations\n')
    reduced <- reduceConfs(merged)
    rm(merged)
    reduced$key <- entry$key
    cat(entry$key, 'Done reducing configurations\n')
    cat(entry$key, 'reduce time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    ##powerStats(reduced$edges)
    p = powerTime(reduced$sched)
    reduced$maxPower = max(p[, power])
    ## this is not entirely precise
    minPower = min(reduced$edges[power > 0, power]) * entry$ranks
    cat(entry$key, 'power stats time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Saving\n')
    save(measurementCols, reduced, entrySpace, countedEntryspace,
         entryCols, entries, confSpace, confCols,
         entry,
         file=filename)
    cat(entry$key, 'Done saving\n')
    cat(entry$key, 'save time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Writing timeslices\n')
    writeSlices(reduced)
    cat(entry$key, 'Done writing timeslices\n')
    cat(entry$key, 'timeslice write time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    ##return(reduced)
  }
  setkeyv(entries, entryCols)
  ##!@todo launch these as separate jobs
  ##result <<-
  rowApply(entrySpace, f)
  ##names(result) <<- entrySpace$key
  NULL
}

f = function(){
  a <<-
    result[[1]]$reduced[is.na(name) & duration > .02,
                                list(hash, OMP_NUM_THREADS, cpuFreq, duration,
                                     pkg_w, pp0_w)][,lapply(.SD, mean),
                                                    by=list(hash, OMP_NUM_THREADS, cpuFreq)][
                                                      order(hash, OMP_NUM_THREADS, cpuFreq)]
  b <<- a[hash==unique(hash)[2]]
  m <<- nnapply(unique(b$cpuFreq),
    function(f){
      b = b[cpuFreq==f,list(OMP_NUM_THREADS, pkg_w,pp0_w)]
      pkg_m = lm(pkg_w ~ OMP_NUM_THREADS, data=b)
      pp0_m = lm(pp0_w ~ OMP_NUM_THREADS, data=b)
      list(pkg=pkg_m,pp0=pp0_m)
    })
}

if(!interactive())
  go()

