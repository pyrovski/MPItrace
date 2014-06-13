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
    result$compEdges[[col]] = entry[[col]]
    if(!any(is.na(result$messageEdges)))
      result$messageEdges[[col]] = entry[[col]]
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
       compEdges=compEdges)
}

## combine within confCols combinations. This will combine multiple
## runs with the same configuration.
reduceConfs = function(x){
  startTime = Sys.time()

  by = c('uid',confCols)
  cat('Reducing between configs\n')

  ## all message edges should be identical between runs
  if(any(is.na(x$messageEdges))){
    cat('No message edges in at least one run\n')
    x$messageEdges = NULL
  } else {
    cat('Merging message edges\n')
    x$messageEdges = rbindlist(x$messageEdges)
    cat('Message edges:', nrow(x$messageEdges), '\n')
    by = c('s_uid','d_uid',confCols)
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
    x$messageEdges[, type:='message']
    x$messageEdges =
      reduceNoEffect(x$messageEdges, c('weight'),
                     setdiff(names(x$messageEdges),
                             c('weight')), c('s_uid','d_uid',confCols))
  }
  
  cat('Computation edges\n')
  x$compEdges = rbindlist(x$compEdges)
  by = c('s_uid',confCols)
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
  x$compEdges[, type:='comp']
  setkey(x$compEdges, d_uid)
  x$compEdges =
    reduceNoEffect(x$compEdges, c('weight','power'),
                   setdiff(names(x$compEdges),
                           c('weight','power')), c('s_uid','d_uid',confCols))

  if(!is.null(x$messageEdges)){
    commonNames = intersect(names(x$messageEdges), names(x$compEdges))
    messageOnlyNames = setdiff(commonNames, names(x$messageEdges))
    compOnlyNames = setdiff(commonNames, names(x$compEdges))
    setkeyv(x$messageEdges, commonNames)
    setkey(x$compEdges, NULL)
    x$edges = merge(x$messageEdges, x$compEdges, all=T)
  } else
    x$edges = x$compEdges
  x$compEdges = NULL
  x$messageEdges = NULL

  ## set power to zero for message edges
  x$edges[is.na(power), power:=0]

  ## assign edge uids
  if(!'e_uid' %in% names(x$edges)){
    setkey(x$edges, s_uid, d_uid, type)
    uids = x$edges[, list(e_uid=.GRP), keyby=list(s_uid, d_uid, type)]
    e = x$edges[uids]
    x$edges = e
    rm(uids, e)
  }
  
  cat('Pareto frontiers\n')
  ## get pareto frontiers
  x$edges = pareto(x$edges)

  ## Get an initial schedule, starting with minimum time per task.
  x$schedule = x$edges[,.SD[which.min(weight)],by=e_uid]

  cat('Schedule and critical path\n')
  schedule = getSchedule(x$schedule)
  ##g = schedule$g
  x$schedule = schedule$edges
  x$critPath = schedule$critPath
  x$vertices = schedule$vertices
  rm(schedule)
  
  ## copy start time to edges
  ## setkey(x$edges, e_uid)
  ## setkey(x$schedule, e_uid)
  ## x$edges = x$edges[x$schedule[, list(e_uid, start)]]

  ## Insert slack edges. These edges will have power, but not minimum
  ## time. To insert the new edges, we need new vertices and new edge
  ## uids. We use the negative of the original edge uid for each.
  cat('Slack edges\n')
  x$schedule = slackEdges(x$schedule, x$critPath)
  x$schedule[is.na(weight), weight:=0]

  ## add slack vertices to vertices table
  x$slackVertices =
    x$schedule[is.na(s_uid),
               list(start=head(start, 1),via=-head(src, 1)), by=src]
  setnames(x$slackVertices, c('vertex', 'start', 'via'))
  ##x$vertices = rbind(x$vertices, slackVertices)
  return(x)
}

writeSlices = function(x, sliceDir='csv'){
  options(scipen=7)
  firstCols = c('e_uid', confCols)
  confName = gsub('[/.]', '_', x$key)
  slices = timeslice(x$schedule, x$edges)
  i = 1
  for(slice in slices){
    sliceName = paste(confName, i, sep='.')
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

    ##!@todo this is easier to get from x$schedule
    write.table(slice[,list(src=head(src, 1), dest=head(dest, 1),
                            edge_rank=head(rank,1)),
                      by=e_uid],
                file=file.path(sliceDir, paste(sliceName, '.edges.csv', sep='')),
                row.names=F, quote=F, sep=',')
    write.table(slice[,c(firstCols, 'weight', 'power'),with=F],
                file=file.path(sliceDir, paste(sliceName, '.edge_weights.csv', sep='')),
                row.names=F, quote=F, sep=',')
    write.table(slice[, list(minPower=min(power), maxPower=max(power)), by=e_uid],
                file=file.path(sliceDir, paste(sliceName, '.edge_powerRange.csv', sep='')),
                row.names=F, quote=F, sep=',')
    write.table(slice[, list(minTime=min(weight), maxTime=max(weight)), by=e_uid],
                file=file.path(sliceDir, paste(sliceName, '.edge_timeRange.csv', sep='')),
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
    i = i + 1
  }
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
    ## startTime = Sys.time()
    ## powerStats(reduced$edges)
    ## cat(entry$key, 'power stats time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Writing timeslices\n')
    writeSlices(reduced)
    cat(entry$key, 'Done writing timeslices\n')
    cat(entry$key, 'timeslice write time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Saving\n')
    save(measurementCols, reduced, entrySpace, countedEntryspace,
         entryCols, entries, confSpace, confCols,
         entry,
         file=filename)
    cat(entry$key, 'Done saving\n')
    cat(entry$key, 'save time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
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

