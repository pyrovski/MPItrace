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
  names(result) = dates

### Merge between runs of the same config. This entails matching
### hashes, requests, and communicators between runs. After the
### matching, the runs should have identical per-rank event
### ordering. UIDs will not match between runs unless we force them to
### match.

  ## Sanity test
  if(length(unique(sapply(result, function(r) nrow(r$compEdges)))) > 1){
    errMsg = 'Runs differ in number of edges!'
    cat(errMsg, '\n')
    edgeTable = sort(table(sapply(result, function(r) nrow(r$compEdges))))
    maxEdgeCount = names(edgeTable[1])
    keep = sapply(result, function(r) nrow(r$compEdges) == maxEdgeCount)
    warning('removing ', length(which(!keep)), ' experiments', immediate.=T)
    result = result[keep]
  }

  f = function(r, name){
    r[[name]]$date = r$date
    r[[name]]
  }
  
  assignments = .rbindlist(lapply(result, f, name='assignments'))
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
  names(collectives) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$collectives=NULL
  gc()

  vertices = lapply(result, '[[', 'vertices')
  names(vertices) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$vertices=NULL

  globals = lapply(result, '[[', 'globals')
  names(globals) = lapply(result, '[[', 'date')
  for(i in 1:length(result)) result[[i]]$globals=NULL
  
  ##!@todo match hashes in vertices tables
  
 
### Comms should already be unified; I mapped MPI_COMM_WORLD and
### MPI_COMM_NULL to -1 and -2, respectively. The other comms should
### align, but we need to test this.

###!@todo I disabled these because we don't have enough RAM to load
###!the runtimes table for all runs in a sequence.

  ## with no created communicators, the comms tables will all be NULL
  
  ## setkey(result, date, rank)
  ## commSeqs = lapply(dates, function(d){
  ##   result[J(d)][comm != '(nil)']$comm
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

  rm(result)

  list(##runtimes=runtimes,
       assignments=assignments,
       messageEdges=messageEdges,
       compEdges=compEdges,
       collectives=collectives,
       vertices=vertices,
       globals=globals)
}

## combine within confCols combinations. This will combine multiple
## runs with the same configuration.
reduceConfs = function(x){
  startTime = Sys.time()

  by = c('uid',confCols)
  cat('Reducing between configs\n')

  x$collectives = base::unique(x$collectives)
  if(length(x$collectives) > 1){
    stop('Unmatched collectives\n')
  }
  x$collectives = x$collectives[[1]]

  x$vertices = base::unique(x$vertices)
  if(length(x$vertices) > 1){
    warning('Unmatched hashes\n', immediate.=T)
    if(length(unique(lapply(x$vertices, function(v) v[, list(vertex, label)]))) > 1){
      stop('Unmatchable hashes\n')
    }
  }
  x$vertices = x$vertices[[1]]

  ##!@todo check for multiple definitions of globals
  x$globals = x$globals[[1]]

  ## all message edges should be identical between runs
  if(any(is.na(x$messageEdges))){
    cat('No message edges in at least one run\n')
    x$messageEdges = NULL
  } else {
    cat('Merging message edges\n')
    measurementCols = c('weight')
### this unique is ok because message weights are identical between
### runs
    x$messageEdges = unique(.rbindlist(x$messageEdges))
    cat('Message edges:', nrow(x$messageEdges), '\n')
    by = c('s_uid','d_uid')
    cores = getOption('mc.cores')
    if(!is.null(cores) && cores > 1){
      setkey(x$messageEdges, s_uid)
      s_uids = unique(x$messageEdges[, s_uid])
      s_uid_chunks = chunk(s_uids, length(s_uids)/cores)
      rm(s_uids)
      x$messageEdges = .rbindlist(mclapply(s_uid_chunks, function(s_uids)
        x$messageEdges[J(s_uids)][,lapply(.SD, mean),by=by]))
      rm(s_uid_chunks)
    } else
      x$messageEdges = x$messageEdges[,lapply(.SD, mean),by=by]
    x$messageEdges_inv =
      x$messageEdges[, setdiff(names(x$messageEdges), measurementCols), with=F]
    x$messageEdges = x$messageEdges[, c('s_uid', 'd_uid', measurementCols), with=F]
  }
  
  cat('Computation edges:', sum(sapply(x$compEdges, nrow)), '\n')

  # extract active thread count from flags bits 4:11
  #
  #!@todo this overwrites OMP_NUM_THREADS from the environment; if we
  #!want this later, we'll have to do something else here
  for(i in 1:length(x$compEdges)){
    x$compEdges[[i]]$OMP_NUM_THREADS =
      pmax(bitwShiftR(bitwAnd(x$compEdges[[i]]$flags, 0xff0), 4), 1)
    x$compEdges[[i]]$flags = bitwAnd(x$compEdges[[i]]$flags, 0x7ffff00f)
  }
  
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
    .rbindlist(mclapply(x$compEdges, function(compEdges)
                       compEdges[, c(by, measurementCols), with=F]
                       ))
  setkey(x$compEdges, s_uid)
  cores = getOption('mc.cores')
  if(!is.null(cores) && cores > 1){
    s_uids = unique(x$compEdges[, s_uid])
    s_uid_chunks = chunk(s_uids, length(s_uids)/cores)
    rm(s_uids)
    x$compEdges = .rbindlist(mclapply(s_uid_chunks, function(s_uids)
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

  ## export an edge set for naive configs

  ##!@todo find a set for naive->thread reduction w/FL.  This is not an
  ##!efficient frontier, but is easily implementable in a runtime
  ##!system.

  ##!@todo consider saving this separately; it is huge relative to the
  ##other edge tables
  x$n_edges = x$compEdges[, .SD[OMP_NUM_THREADS==max(OMP_NUM_THREADS)], by=e_uid]

  ## plot perf vs power for one edge
  setkey(x$compEdges, e_uid)
  setkey(x$edges_inv, e_uid)
  maxThreads = max(x$compEdges$OMP_NUM_THREADS, na.rm=T)
  threadEdges = unique(x$compEdges[OMP_NUM_THREADS == maxThreads, e_uid])
  minWeightEdges = x$compEdges[J(threadEdges)][x$edges_inv[,list(e_uid, rank)]][,
    list(minWeight=min(weight), rank),by=e_uid][minWeight > .02][,
                                        head(.SD, 4),by=rank][, e_uid]
  plotEdges = minWeightEdges
  if(length(plotEdges)){
    for(e in plotEdges)
      plotPerfPower(x$compEdges[J(e)],
                    name=paste(x$key, 'Phase', e))
  }
  
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
  x$vertices$newVertex = 1:nrow(x$vertices)
  setkey(x$vertices, vertex)
  x$edges_inv = data.table::copy(x$edges_inv) ## must have been a faulty rbindlist somewhere
  for(col in c('src', 'dest', 'o_dest')){
    setkeyv(x$edges_inv, col)
    toMatch = unique(Filter(Negate(is.na), x$edges_inv[[col]]))
    x$edges_inv[x$vertices[J(toMatch), list(vertex, newVertex)], c(col) := list(newVertex)]
  }
  x$vertices[, c('vertex', 'newVertex') := list(newVertex, NULL)]
  
  ## Get an initial schedule, starting with minimum time per task.
  x$schedule = x$edges[,.SD[which.min(weight)],keyby=e_uid]
  ## get src and dest vertices
  x$schedule = x$schedule[x$edges_inv[, list(e_uid, src, dest, rank, type)]]
  
  cat('Schedule and critical path\n')
  schedule = getSchedule(x$schedule, x$vertices)
  x$schedule = schedule$edges
  x$critPath = schedule$critPath
  x$vertices = schedule$vertices
  x$minTime = x$vertices[vertex==2,start]
  needSlack = schedule$needSlack
  rm(schedule); gc()

  #! relax schedule
  #!@todo the join with vertices could be done for all edges at once
  setkey(x$edges, e_uid)
  f = function(a, b) a <= b + 1e-8
  x$schedule = x$schedule[,
    {
      deadline=x$vertices[J(dest), start];
      wslack=deadline-start;
      e=x$edges[J(.BY[[1]])][f(weight, wslack)][which.max(weight)];
      if(!nrow(e)){
        stop('expected at least one row')
      }
      a=.SD[,setdiff(names(.SD),names(e)),with=F];
      e=cbind(a,e);
      e[, e_uid:=NULL]
    },
  by=e_uid]
  
  ## Insert slack edges. These edges will have power, but not minimum
  ## time. To insert the new edges, we need new vertices and new edge
  ## uids. We use the negative of the original edge uid for each.
  cat('Slack edges\n')
  startTime = Sys.time()
  x$schedule = x$schedule[x$edges_inv[, list(e_uid, s_uid, d_uid)]]
  activeWaitConf = x$edges[power > 0][which.min(power), c(confCols, 'power'), with=F]
  activeWaitConf$weight = as.numeric(NA)
  x$schedule = slackEdges(x$schedule, activeWaitConf, x$critPath, needSlack)
  x$schedule[is.na(weight), weight:=0]

  ## add slack vertices to vertices table
  slackVerticesPost =
    x$schedule[is.na(s_uid) & !is.na(d_uid), ## post
               list(start=head(start, 1),via=-head(src, 1)),
               by=src]
  setnames(slackVerticesPost, c('vertex', 'start', 'via'))

  if(attr(x$schedule, 'doubleSlack')){
    slackVerticesPre =
      x$schedule[!is.na(s_uid) & is.na(d_uid), ## pre
                 list(start=head(start, 1),via=as.numeric(NA)),
                 by=dest]
    setnames(slackVerticesPre, c('vertex', 'start', 'via'))
    x$slackVertices = rbind(slackVerticesPre, slackVerticesPost)
    rm(slackVerticesPost, slackVerticesPre)
  } else
    x$slackVertices = slackVerticesPost
  rm(slackVerticesPost)
  
  #!@todo fix?
  x$slackVertices[ ,c('hash', 'label', 'reqs') := as.character(NA)]
  
  setcolorder(x$slackVertices, names(x$vertices))
  cat(x$key, 'Slack time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')

  writeGraphFromSchedule(
    x$schedule, x$critPath, x$key, compile=F,
    v=rbind(x$vertices,x$slackVertices)[,
      list(vertex,label)][,label:=paste(label,vertex,sep=': ')])
  return(x)
}

writeSlices = function(x, sliceDir='csv'){
  dir.create(sliceDir, showWarnings=F)
  options(scipen=7)
  firstCols = c('e_uid', confCols)
  confName = gsub('[/. ]', '_', x$key)
  files = Sys.glob(paste(file.path(sliceDir,confName), '*.csv', sep=''))
  mclapply(files, unlink)
  
  schedVertices = rbind(x$vertices, x$slackVertices)
###!@todo this would use less memory if we computed one slice at a time
  ##slices = timeslice(x$schedule, schedVertices, x$edges)
  ##names(slices) = sprintf('%.3f', as.numeric(names(slices)))

  writeSlice = function(slice, sliceTime, schedule, schedVertices, reachable){
    setkey(slice, e_uid)
    sliceName = paste(confName, sliceTime, sep='_')
    setcolorder(slice, c(firstCols, setdiff(names(slice), firstCols)))

    writeGraphFromSchedule(
      schedule, key=sliceName, compile=F,
      v=schedVertices[, label:=paste(label,vertex,sep=': ')])
    
    if(length(grep('ILP', sliceTime)))
      sliceNameFixed = paste(confName, sub('ILP', 'fixedLP', sliceTime), sep='_')
    
    writeTable = function(table, suffix){
      base = paste(sliceName, suffix, sep='')
      filename = file.path(sliceDir, base)
      write.table(table, file=filename,
                  quote=F, sep=',', row.names=F)
      if(length(grep('ILP', sliceTime))){
        filenameFixed = file.path(sliceDir, paste(sliceNameFixed, suffix, sep=''))
        if(!file.exists(filenameFixed))
           file.symlink(base, filenameFixed)
      }
    }

    ## write confSpace
    writeTable(unique(slice[,confCols,with=F], by=confCols), '.confSpace.csv')

    if(!length(grep('ILP', sliceTime))){
      writeTable(slice[,list(vertex=union(src,dest))], '.vertices.csv')
      
      tmpSlice =
        slice[,list(
          ##!@todo this is easier to get from schedule
          src=head(src, 1), dest=head(dest, 1),
          edge_rank=head(rank,1),
          ## but not these
          minTime=min(weight),
          maxTime=max(weight),
          minPower=min(power),
          maxPower=max(power),
          left=head(left,1),
          right=head(right,1)),
              by=e_uid]
    } else { ## ILP
      base = paste(sliceName, '.reachable.dat', sep='')
      filename = file.path(sliceDir, base)
      filenameFixed =
        file.path(sliceDir, paste(sliceNameFixed, '.reachable.dat', sep=''))
      reachableFile = file(filename, 'w')
      for(vertex in names(reachable))
        write(paste(vertex, ': ', paste(reachable[[vertex]], collapse=' '),
                    sep=''), reachableFile)
      close(reachableFile)
      if(!file.exists(filenameFixed))
        file.symlink(base, filenameFixed)
      
      writeTable(schedVertices[, list(vertex, ancestors, descendants,
                                      vertexEvent=vertexOrder)], '.vertices.csv')

      ## write precedence matrix
      e2 = schedule[,list(e_uid, src)]
      setkey(e2, src)
      setkey(schedule, e_uid)
      precedence =
        .rbindlist(lapply(schedule[, e_uid],
                          function(e){
                            d = schedule[J(e), dest]
                            succ=e2[J(d), list(successor=e_uid)]
                            succ$edge=e;succ
                          }))
      setcolorder(precedence, c('edge', 'successor'))
      writeTable(precedence, '.precedence.csv')
      
      tmpSlice =
        slice[,list(
          ##!@todo this is easier to get from schedule
          src=head(src, 1), dest=head(dest, 1),
          edge_rank=head(rank,1),
          ## but not these
          minTime=min(weight),
          maxTime=max(weight),
          minPower=min(power),
          maxPower=max(power)),
              by=e_uid]
    }
    writeTable(tmpSlice, '.edges.csv')
    rm(tmpSlice)
    
    writeTable(slice[,c(firstCols, 'weight', 'power'),with=F],
               '.edge_weights.csv')
    writeTable(slice[,list(count=nrow(.SD)), by=e_uid][count > 1, list(e_uid)],
                 '.edge_multiConf.csv')

    ## ##!I only need a single rank column. Slack edges always go on the
    ## ##!destination rank, and we don't care about messages. Only the
    ## ##!computation and slack edges should have rank information.

    ranks = unique(slice[, list(rank)])
    writeTable(ranks, '.rank.csv')
    base = paste(sliceName, '.Rsave', sep='')
    filename = file.path(sliceDir, base)
    save(slice, file=filename)
    if(length(grep('ILP', sliceTime))){
      filenameFixed =
        file.path(sliceDir, paste(sliceNameFixed, '.Rsave', sep=''))
      if(!file.exists(filenameFixed))
         file.symlink(base, filenameFixed)
    }
  }
  ## after the fixed-order LP, we don't need this anymore.
  ##mclapply(names(slices), function(sliceTime)
  ##         writeSlice(slices[[sliceTime]], sliceTime, x$schedule))
  setkey(x$edges, e_uid)
  setkey(x$schedule, e_uid)
  
  ## add slack edges to table passed to writeSlice for ILP
  if(any(x$schedule[, e_uid] < 0)){
    result =
      rbind(x$edges[x$schedule[e_uid > 0, list(e_uid)]],
            x$schedule[e_uid < 0, names(x$edges), with=F])
  } else
    result = x$edges[x$schedule]
  setkey(result, e_uid)
  result = result[x$schedule[, list(e_uid, src, dest, rank, type)]]

  setkey(schedVertices, vertex)
  ## ancestors, descendants
  g = graph.data.frame(x$schedule[, list(src, dest)])
  vertexNames = V(g)$name
### Ancestors of each vertex
  ancestors = neighborhood.size(g, order=vcount(g), mode='in') - 1
  schedVertices$ancestors = ancestors[order(as.numeric(vertexNames))]
### Descendants of each vertex
  descendants = neighborhood.size(g, order=vcount(g), mode='out') - 1
  schedVertices$descendants = descendants[order(as.numeric(vertexNames))]
  rm(ancestors, descendants)
  reachable = neighborhood(g, order=vcount(g), mode='out')
  names(reachable) = vertexNames
  reachable =
    napply(reachable, function(x, name)
           Filter(function(x) x != name, vertexNames[x]))
  rm(vertexNames)
  
  graphFile = file.path(sliceDir, paste(confName, '.graph.dot', sep=''))
  write.graph(g, file=graphFile, format='dot')
  system(paste('gzip -f ', graphFile, sep=''), wait=F)

  ## compute vertex order
  vo =
    schedVertices[vertex > 0][order(start)][,list(
      vertex, vertexOrder=.GRP-1
      ),by=start][,list(vertexOrder), keyby=vertex]
  schedVertices[vo, vertexOrder:=vertexOrder]
  schedVertices[is.na(vertexOrder), vertexOrder:=-1]
  rm(vo)

### write barrier-separated sections separately. This
### corresponds to finding vertices through which all paths pass. We
### should be able to find such vertices in O(|V|^2*|E|) time by
### independently deleting each vertex and testing the number of
### connected components. We can reduce the list of candidate vertices
### by looking for vertices with outgoing edges to every rank, and
### vertices matching MPI_Barrier, although these qualities are not
### sufficient. The resulting complexity is O(|B|*|V|*|E|), where B is
### the set of barrier/high degree edges.

### don't forget that the first vertex in each section is 1
### and last vertex is 2!

  x$numRanks = length(unique(x$edges_inv$rank))
  cuts = list()
  if(x$numRanks > 1){
    candidates = degree(g, mode='in') >= x$numRanks &
      degree(g, mode='out') >= x$numRanks
    if(any(candidates)){
      candidates = names(which(candidates))
      cuts = lapply(candidates, function(v){
        result = NULL
        if(no.clusters(delete.vertices(g, v)) > 1)
          v
        else
          NA
      })
      cuts = cuts[!is.na(cuts)]
    }
  }
  if(length(cuts)){
### perform the cuts: renumber first/last vertices, write separate file sets
    
### produce a list of data tables of length (number of cuts + 1)
    setkey(schedVertices, vertex)
    # get cuts in reverse order of vertex start time
    cuts = schedVertices[J(as.numeric(cuts))][order(-start)]$vertex
    startVertex = '1'
    endVertex = '2'
    schedule = data.table::copy(x$schedule)
    for(v in cuts){
      ## cut graph at v (remove outgoing edges)
      cutEdges = incident(g, as.character(v), mode='out')
      outNeighbors = neighbors(g, as.character(v), mode='out')
      g2 = delete.edges(g, cutEdges)
      frontVertices = V(g2)$name[subcomponent(g2, startVertex)]
      backVertices = V(g2)$name[subcomponent(g2, endVertex)]

      ## retain only front component
      g2 = induced.subgraph(g, c(v, backVertices))
      ## g2 = g2 + vertex(as.character(v))
      outNeighbors = V(g)$name[outNeighbors]
      g2[from=as.character(rep(v, times=length(outNeighbors))),
         to=outNeighbors] = T
      g = induced.subgraph(g, frontVertices)
      
      ## rename new sink vertex in g to '2'
      V(g)[as.character(v)]$name = '2'

      frontVertices = as.numeric(frontVertices)
      frontVertices = setdiff(frontVertices, v)
      backVertices = as.numeric(backVertices)

      ## rename new source vertex in g2 to '1'
      ##V(g2)[as.character(v)]$name = '1'

### Ancestors of each vertex
      ancestors = neighborhood.size(g2, order=vcount(g2), mode='in') - 1
### Descendants of each vertex
      descendants = neighborhood.size(g2, order=vcount(g2), mode='out') - 1
      backSchedVertices = schedVertices[J(c(v, backVertices))]
      setkey(backSchedVertices, vertex)
      vertexNames = V(g2)$name
      backSchedVertices$ancestors = ancestors[order(as.numeric(vertexNames))]
      backSchedVertices$descendants = descendants[order(as.numeric(vertexNames))]
      rm(ancestors, descendants)
      
      reachable = neighborhood(g2, order=vcount(g2), mode='out')
      names(reachable) = vertexNames
      reachable =
        napply(reachable, function(x, name)
               Filter(function(x) x != name, vertexNames[x]))
      rm(vertexNames)
      
      rm(g2)
      schedVertices = schedVertices[J(c(frontVertices, v))]

      ## filter tables
      backResult = result[src %in% backVertices | dest %in% backVertices]
      result = result[src %in% frontVertices | dest %in% frontVertices]
      backSchedule = schedule[src %in% backVertices | dest %in% backVertices]
      schedule = schedule[src %in% frontVertices | dest %in% frontVertices]

      ## adjust start times
      backSchedule$start = backSchedule$start - min(backSchedule$start)
      
###renumber vertices in backSchedule, backResult, schedule, result
      backSchedule[src==v, src := 1]
      backResult[src==v, src := 1]
      backSchedVertices[vertex==v, vertex := 1]
      schedVertices[vertex==v, vertex := 2]
      setkey(schedVertices, vertex)
      schedule[dest==v, dest := 2]
      result[dest==v, dest := 2]
      names(reachable)[names(reachable) == as.character(v)] = '1'
      
      writeSlice(backResult, paste('ILP.cut_', v, sep=''), backSchedule,
                 backSchedVertices, reachable=reachable)
    }
    ## write earliest chunk
    ## Ancestors of each vertex
    ancestors = neighborhood.size(g, order=vcount(g), mode='in') - 1
    ## Descendants of each vertex
    descendants = neighborhood.size(g, order=vcount(g), mode='out') - 1
    setkey(schedVertices, vertex)
    vertexNames = V(g)$name
    schedVertices$ancestors = ancestors[order(as.numeric(vertexNames))]
    schedVertices$descendants = descendants[order(as.numeric(vertexNames))]
    rm(ancestors, descendants)
    reachable = neighborhood(g, order=vcount(g), mode='out')
    names(reachable) = vertexNames
    reachable =
      napply(reachable, function(x, name)
             Filter(function(x) x != name, vertexNames[x]))
    rm(g)

    writeSlice(result, 'ILP.cut_1', schedule, schedVertices,
               reachable=reachable)
  } else ## no cuts
    writeSlice(result, sliceTime = 'ILP.cut_1', x$schedule, schedVertices,
               reachable=reachable)
  confName
}

##!@todo this may run into memory limitations. If so, just run the
##!configurations one at a time.

go = function(force=F){
  ## group experiments by entryCols from entries

  load('mergedEntries.Rsave', envir=.GlobalEnv)
  ##! remove runs with turboboost
  entries = entries[cpuFreq != 2601000]
  if(!nrow(entries)){
    cat('no entries!\n')
    return()
  }
  setkeyv(entries, entryCols)
  
  measurementCols <<- c('duration','pkg_w','pp0_w','dram_w')
  confSpace <<- unique(entries[,confCols,with=F])
  setkey(confSpace)

  f = function(entry){
    startTime = Sys.time()
    filename = paste('mergedData', gsub('[/. ]', '_', entry$key),
      'Rsave', sep='.')    
    if(file.exists(filename) & !force){
      cat(entry$key, 'already merged\n')
      return(NULL)
    }

### check for missing configs
    confs = unique(entries[entry, confCols, with=F])
    if(nrow(confs) == 0){
      warning('no entries for ', entry, '\n', immediate.=T)
      return(NULL)
    }
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
    merged$key = entry$key
    reduced <- reduceConfs(merged)
    rm(merged)
    cat(entry$key, 'Done reducing configurations\n')
    cat(entry$key, 'reduce time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    ##powerStats(reduced$edges)
    ##p = powerTime(reduced$sched, rbind(reduced$vertices, reduced$slackVertices))

    setkey(reduced$edges, e_uid)
    setkey(reduced$edges_inv, e_uid)
    reduced$maxPower =
      sum(reduced$edges[,list(e_uid, power)][,.SD[which.max(power)],
                                             by=e_uid][reduced$edges_inv[,
                                               list(e_uid, rank)]][,
                                                                   .SD[which.max(power)],
                                                                   by=rank]$power)

    ## this is non-slack power; slack edges are not present in reduced$edges.
    reduced$minPower =
      sum(reduced$edges[,list(e_uid, power)][,.SD[which.min(power)],
                                             by=e_uid][reduced$edges_inv[,
                                               list(e_uid, rank)]][,
                                                                   .SD[which.max(power)],
                                                                   by=rank]$power)
    cat(entry$key, 'power stats time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Saving\n')
    save(measurementCols, reduced, entrySpace, countedEntryspace,
         entryCols, entries, confSpace, confCols,
         entry, missingConfs,
         file=filename)
    cat(entry$key, 'Done saving\n')
    cat(entry$key, 'save time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    startTime = Sys.time()
    cat(entry$key, 'Writing timeslices\n')
    writeSlices(reduced)
    cat(entry$key, 'Done writing timeslices\n')
    cat(entry$key, 'timeslice write time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
    ##return(reduced)
    NULL
  }
  setkeyv(entries, entryCols)
  ##!@todo launch these as separate jobs
  rowApply(entrySpace, f)
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

