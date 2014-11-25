##!@todo finish
mcdtby = function(x, chunkBy, by, f, SDcols='all'){
  cores = getOption('mc.cores')
  if(length(chunkBy) != 1)
    stop('expected chunkBy of length 1\n')
  if(!is.null(cores) && cores > 1){
    setkeyv(x, by)
    groups = unique(x[, chunkBy, with=F])
    sets = chunk(groups, length(groups)/cores)
    ##!@todo
    rm(groups)
    result = .rbindlist(mclapply(sets, function(s) x[s]))
  } else {
    if(SDcols == 'all')
      result = x[, f(.SD), by=by]
    else
      result = x[, f(.SD), by=by, .SDcols=SDcols]
  }
  return(result)
}

oldPowerTime = function(x){
  ranks  = unique(x$rank)
  powerCols = grep('_w$',names(x),value=T)
  cols = c('start','duration',powerCols)
  powerSteps =
    lapply(ranks, function(r){
      x = x[rank == r][order(uid), cols, with=F]
      x = x[,list(start, power=pkg_w+dram_w)] ## pp0_w is included in pkg_w
      list(times = x$start, steps = stepfun(x$start, c(x$power,tail(x$power,1))))
    })
  times = sort(unique(unlist(lapply(powerSteps, '[[', 'times'))))
  steps = lapply(powerSteps, '[[', 'steps')
  powers =
    as.data.table(cbind(start=times,
          power=rowSums(as.data.table(lapply(steps, function(s) s(times))))))
  return(powers)
}

powerTime = function(edges, vertices){
  ranks = unique(edges[, rank])
  f = function(x) c(x, 0)
  rf = function(r){
    edges = edges[rank == r & power > 0]
    ##! data.table doesn't order correctly?
    ##[order(start)]
    o = order(edges[,start])
    edges = edges[o]
    maxStart = max(vertices[, start])
    times = edges[, start]
    
    powers = data.table(start=times, power=edges[, power])
    cs = cumsum(rle(rev(powers[, power]))$lengths)
    sel = rev(tail(cs, 1) + 1 - cs)
    powers = powers[sel]
    
    steps = stepfun(x=c(powers$start, maxStart), y=c(0, powers$power, 0))
    list(times=times,
         ##powers=powers,
         steps=steps)
  }
  powerSteps = mclapply(ranks, rf)
  times = c(-.000001, sort(unique(unlist(lapply(powerSteps, '[[', 'times')))))
  steps = lapply(powerSteps, '[[', 'steps')
  powers =
    as.data.table(cbind(start=times,
          power=rowSums(as.data.table(lapply(steps, function(s) s(times))))))

  ## remove runs of identical powers
  cs = cumsum(rle(rev(powers[, power]))$lengths)
  sel = rev(tail(cs, 1) + 1 - cs)
  powers = powers[sel]
  
  return(powers)
}

plotPowerTime = function(pt, name){
  pdf(file=paste(name, '.pdf', sep=''))
  s = stepfun(pt[[1]], c(pt[[2]], 0))
  plot(s, main=name)
  dev.off()
  NULL
}

rebuildScheduleWithSlack = function(e, activeWaitConf){
  sched = getSchedule(e)
  schedVertices = sched$vertices
  sched = slackEdges(sched$edges, activeWaitConf, sched$critPath)
  slackVertices =
    sched[is.na(s_uid),
          list(start=head(start, 1),via=-head(src, 1)), by=src]
  setnames(slackVertices, c('vertex', 'start', 'via'))

  ## nullify 0-length slack edges
  ranks = unique(sched[, rank])
  setkey(sched, rank)
  sched = .rbindlist(mclapply(ranks, function(r){
    sched = sched[J(r)][order(start, -e_uid)]
    dupes = duplicated(sched[, start])
    sched[dupes, power := 0]
    sched
  }))

  vertices = rbind(schedVertices, slackVertices)
  pt = powerTime(sched, vertices)
  list(pt=pt, sched=sched, vertices=vertices)
}

powerStats = function(edges, edges_inv, powerLimits, limitedOnly=F, name){
### need to assign start times for each set 
  ##!@todo measure this in runtime in post-init
  ##!@todo this should be processor-specific
  activeWaitConf = edges[power > 0][which.min(power), c(confCols, 'power'), with=F]
  activeWaitConf$weight = as.numeric(NA)
  
  setkey(edges, e_uid)
  setkey(edges_inv, e_uid)
  
  cols = c('e_uid', 's_uid', 'd_uid', 'src', 'dest', 'rank', 'type')
  cols = intersect(cols, names(edges_inv))

  f = function(e){
    e = e[edges_inv[, cols, with=F]]
    result = rebuildScheduleWithSlack(e, activeWaitConf)
    result$pt
  }

  if(!limitedOnly){
    ##fastest
    powersMinTime =
      mcparallel(f(edges[,.SD[which.min(weight)],by=e_uid]), name='powersMinTime')
    
    ## max power
    powersMaxPower =
      mcparallel(f(edges[,.SD[which.max(power)],by=e_uid]), name='powersMaxPower')

    ## min power
    powersMinPower =
      mcparallel(f(edges[,.SD[which.min(power)],by=e_uid]), name='powersMinPower')

    mcresult = mccollect()
    
    result =
      list(powersMinTime = mcresult[[as.character(powersMinTime$pid)]],
           powersMaxPower = mcresult[[as.character(powersMaxPower$pid)]],
           powersMinPower = mcresult[[as.character(powersMinPower$pid)]])
    
    ## result =
    ##   list(minTime  = powersMinTime,
    ##        maxPower = powersMaxPower,
    ##        minPower = powersMinPower)
  } else
    result = list()
  
  ## power-limited
  if(!missing(powerLimits)){
    ranks = unique(edges_inv[, rank])
    nRanks = length(ranks)
###!@todo implement a few strategies for runtime power management: 1:
### uniform per-rank power constraints, naive config+FL 2: uniform
### per-rank power constraints, intelligent config 3: initially
### uniform per-rank power constraints with excess reallocation (this
### requires event simulation).
    result_pl =
      mclapply(powerLimits, function(powerLimit){
        ##!@todo get maxconf.  This could depend on all confCols.
        
        ## this may violate the power constraint. !@todo warn on violation
        
        ##!@todo retain inefficient configurations for pl1
        
        pl2 =
          edges[, .SD[power <= powerLimit/nRanks |
                      power == min(power)][which.min(weight)],
                by=e_uid]
        
        pl2 = f(pl2)
      })
    names(result_pl) = paste(powerLimits, 'pl2', sep='_')
    result = c(result, result_pl)
    rm(result_pl)
  }
  if(!missing(name)){
    names(result) = paste(name, names(result))
  }
  napply(result, plotPowerTime, mc=T)
  return(result)
}

## return sets of indices for each timeslice and how much of each task
## is in each timeslice. Slack edges have zero weight, so they will
## never be bisected by a timeslice. They are important to include
## because they contribute to power consumption.

timeslice = function(sched, vertices, edges,
  ts_start=0, length=.01, n=Inf)
{
  setkey(vertices, vertex)
###! this only looks at the start time of the next edge because we
###! force timeslices to match the original schedule
  setkey(sched, dest)
  sched = sched[vertices[, list(vertex, deadline=start)]]
  if(nrow(sched[start > deadline]) > 0)
    stop('invalid deadline!\n')
  ##!@todo start time will change if we schedule from the back
  sched[, wslack:=deadline-start]
  maxDeadline = max(sched[,deadline])
  if(ts_start > maxDeadline)
    stop('Start past last task')

  ranks = unique(sched[, rank])
  
  sliceTimes = head(seq(ts_start, maxDeadline, length), n)
  setkey(sched, e_uid)
  setkey(edges, e_uid)

  f = function(sliceTime){
    nextSlice = sliceTime + length
### There's a package for finding interval overlap, but we only do
### this once.

    ## interior edges: entirely included in the interval
    int =
      sched[start >= sliceTime & deadline < nextSlice,
            list(e_uid, weight, wslack, right=wslack)]
    int[, c('left', 'splitSrc', 'splitDest') := list(0, F, F)]

    ## exterior edges: edges include interval
    ##! re-number src and dest vertices
    ext =
      sched[start <  sliceTime & deadline >= nextSlice,
            list(e_uid, weight, wslack,
                 left=sliceTime-start,
                 right=nextSlice-start)]
    ext[, c('splitSrc', 'splitDest') := list(T, T)]

    ##! re-number src vertex
    ## left overlap
    left =
      sched[start < sliceTime & deadline > sliceTime & deadline < nextSlice,
            list(e_uid, weight, wslack,
                 left=sliceTime-start, right=wslack)]
    left[, c('splitSrc', 'splitDest') := list(T, F)]
    
    ## right overlap
    ##! re-number dest vertex
    right =
      sched[start >= sliceTime & start < nextSlice & deadline >= nextSlice,
            list(e_uid, weight, wslack,
                 right=nextSlice-start)]
    right[, left := 0]
    right[, c('splitSrc', 'splitDest') := list(F, T)]
    
    result = rbind(int, ext, left, right, use.names=T)
    result[weight == 0, c('left', 'right') := list(0,0)]
    ## at this point, right and left are in seconds. We want them to
    ## be fractions of the original edge weight.
    if(nrow(result[left > right]) > 0){
      stop('invalid slicing\n')
    }

###!we need to modify all edges according to the fraction in the
###!timeslice, not just the scheduled edges. So, find out which edges
###!are in the timeslice, then determine their splits. Possible splits
###!depend on the type of overlap.

    ##,'left','right'
    ##, NULL, NULL
    result[weight != 0, c('frac', 'left', 'right') :=
           list((right-left)/wslack, left/wslack, right/wslack)]
    result[weight == 0, frac := 0]
    ##, 'wslack'
    result[, c('weight') := list(NULL)]
    setkey(result, e_uid)
    if(any(result[, e_uid] < 0)){
      ## get multiple configs for non-slack edges, one config for slack edges
      result =
        rbind(edges[result[e_uid > 0]],
              sched[, names(edges), with=F][result[e_uid < 0]])
    } else
      result = edges[result]
    result[,c('weight', 'oWeight', 'frac') := list(frac*weight, weight, NULL)]
    if(nrow(result[weight < 0]) > 0)
      stop('Negative-weight edge(s)\n')
    setkey(result, e_uid)
    result = result[sched[, list(e_uid, src, dest, rank, type)]]
    if(any(!ranks %in% unique(result[, rank])))
      stop('every rank should have an edge in every timeslice\n')
    return(result)
  }
  
  result = nnapply(sliceTimes, f, mc=T)
  result
}

minConf = function(confs){
  ## if we can't simultaneously minimize all knobs, choose the config
  ## with the lowest power. This may be arbitrary for modeled power.
  confs[power > 0][which.min(power)]
}

## orders vertices and edges in topological order, adds a start time
## and e_uid to each edge, and returns the critical path.
getSchedule = function(edges, vertices=edges[,list(vertex=union(src,dest))],
  doCritPath=T){
### this function is intended to be called with one edge row per e_uid
  if(any(edges[, list(count=nrow(.SD)), by=e_uid]$count > 1))
    stop('Duplicate e_uids in getSchedule\n')

  numericVertices = T
  if(inherits(vertices$vertex, 'character'))
    numericVertices = F
  else if(any(diff(sort(vertices[, vertex])) > 1) || !1 %in% vertices[, vertex])
    stop('Vertices must be numbered 1:nrow(vertices)\n')
  
  cat('Graph construction\n')
  setcolorder(edges,
              c('src','dest',setdiff(names(edges), c('src','dest'))))
  if(any(edges[, weight] < 0)){
    stop('Negative-weight edge(s)!')
  }
  g = graph.data.frame(edges[,c('src','dest'), with=F])
  if(no.clusters(g) > 1)
    stop('Graph has more than one cluster!')

  ## this assumes we're going to start edges as soon as their
  ## dependencies are met
  d = degree(g, mode='in') > 1
  needSlack = as.numeric(names(d[which(d)]))
  rm(d)
  
  gd = lapply(get.data.frame(g, what='both'), as.data.table)
  if(numericVertices)
    gd$vertices$name = as.numeric(gd$vertices$name)
  cat('Topological ordering\n')
  ts_order = topological.sort(g)

  ## cat('Vertex ancestry (for ILP taskEvent fixing)\n')
### Ancestors of each vertex
  ## ancestors = neighborhood.size(g, order=vcount(g), mode='in') - 1
  ## ancestors = ancestors[order(gd$vertices$name)]
    
### Descendants of each vertex
  ## descendants = neighborhood.size(g, order=vcount(g), mode='out') - 1
  ## descendants = descendants[order(gd$vertices$name)]

  rm(g)
  
  setkey(vertices, vertex)
  ## vertices$ancestors = ancestors
  ## vertices$descendants = descendants
  ## rm(ancestors, descendants)
  vertices_TO = data.table::copy(vertices[J(gd$vertices[ts_order]), list(vertex)])
  rm(gd); gc()

  if(!numericVertices){
    ## renumber vertices temporarily
    vertexMap = vertices_TO[, list(vertex)]
    vertexMap[, newVertex := .GRP, by=vertex]
    setkey(vertexMap, vertex)
    vertices = vertices[vertexMap]
    vertices$vertex = NULL
    setnames(vertices, 'newVertex', 'vertex')
    setkey(vertices, vertex)

    setkey(vertices_TO, vertex)
    vertices_TO = vertices_TO[vertexMap]
    vertices_TO$vertex = NULL
    setnames(vertices_TO, 'newVertex', 'vertex')
    setkey(vertices_TO, vertex)
    
    setkey(edges, src)
    edges = vertexMap[edges]
    edges$vertex = NULL
    setnames(edges, 'newVertex', 'src')

    setkey(edges, dest)
    edges = vertexMap[edges]
    edges$vertex = NULL
    setnames(edges, 'newVertex', 'dest')
  }

  setkey(edges, src)

  ## define a start time for each edge
  cat('Start times\n')
  startTime = Sys.time()

  vertices$start = -Inf
  vertices[J(vertices_TO[1, list(vertex)]), start := 0]
  if(numericVertices){
    vertices[, via:=as.integer(NA)]
  } else
    vertices[, via:=as.character(NA)]

  count = 0
  progress = -1
  
###!@todo this could be faster if we selected only the relevant set of
###!vertices for each src, then merged after the loop

###!@todo using the igraph query mechanisms may also improve the speed
  ## I checked this, and it seems that the '[[' operator time is the
  ## same regardless of whether I want info for one vertex or many,
  ## and g[[1]] also slower than evaluating e[src==1, dest]. However,
  ## incident() and neighbors() seem to be fast.  neighbors() seems to
  ## take edge direction into account.
  
###!@todo using an LP solver for this would be way faster
  for(vertex in vertices_TO$vertex){
    ## get start time for src vertex
    v = vertices[vertex, list(src=vertex, start)]
    setkey(v, src)      
    outEdges = edges[v, list(src, dest, e_uid, start, weight)]
    if(nrow(outEdges) < 1)
      next
    setkey(outEdges, src)

    ## get potential start times for dest vertices
    outEdges[, c('start', 'start_u') := list(NULL, start+weight)]

    ## resolve cases with multiple edges between src and dest
    outEdges = outEdges[, .SD[which.max(start_u)], keyby=dest]
    ## should now have one edge per dest in v (already have a single src)

    ## update start times for dest vertices
    ##outEdges = vertices[outEdges[, list(dest, start_u, e_uid)]]
    outEdges[, start := vertices[outEdges[,dest], list(start)]]
    ## join on the dest vertex from v and 'vertex' from vertices
    outEdges = outEdges[start_u > start, list(dest, start_u, e_uid)]
    setnames(outEdges, c('vertex', 'start', 'via'))
    vertices[outEdges[,vertex], c('vertex', 'start', 'via') := outEdges]
    
    if(!count %% 1000){
      if(progress != (n_progress <- round(count/nrow(vertices)*100)))
        cat(n_progress, '%\n')
      progress = n_progress
    }
    count = count + 1
  }
  cat(100, '%\n')

  if(!numericVertices){
    ## undo vertex renumbering
    setkey(vertexMap, newVertex)
    
    setnames(vertices, 'vertex', 'newVertex')
    setkey(vertices, newVertex)
    vertices = vertices[vertexMap]
    vertices$newVertex = NULL

    setnames(vertices_TO, 'vertex', 'newVertex')
    setkey(vertices_TO, newVertex)
    vertices_TO = vertices_TO[vertexMap]
    vertices_TO$newVertex = NULL

    setnames(edges, 'src', 'newVertex')
    setkey(edges, newVertex)
    edges = vertexMap[edges]
    edges$newVertex = NULL
    setnames(edges, 'vertex', 'src')
    
    setnames(edges, 'dest', 'newVertex')
    setkey(edges, newVertex)
    edges = vertexMap[edges]
    edges$newVertex = NULL
    setnames(edges, 'vertex', 'dest')

    rm(vertexMap)
  }

###!@todo this assumes edges start as soon as their dependencies are
###!satisfied. If we want to schedule from the back, edges should
###!start as late as possible; for each edge:
###! start = vertices[J(dest), start] - weight
  setkey(vertices, vertex)
##  edges = edges[vertices[, list(vertex, start)]]
  setkey(edges, src)
  edges = vertices[, list(vertex, start)][edges]
  setnames(edges, 'vertex', 'src')
 
  ##if(any(edges$start == -Inf))
  ##  stop('Some edges not assigned a start time!\n')
  cat('Start times time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')

  if(doCritPath){
    cat('Critical path\n')

###!igraph doesn't support finding the shortest path with negative
###!edge weights.

### Does it support finding the number of edges in the shortest path?
### Yes.

    setkey(edges, e_uid)
    setkey(vertices, vertex)
    c_vertex = tail(vertices_TO[, vertex], 1)
    head_vertex = vertices_TO[1, vertex]
    critPath = c()
    while(c_vertex != head_vertex){
      v = vertices[J(c_vertex)][, via]
      c_vertex = edges[J(v), src, mult='first']
      critPath = c(v, critPath)
    }
  } else
    critPath=NULL

  return(list(edges=edges, vertices=vertices,
              critPath=unname(unlist(critPath)),
              needSlack=needSlack))
}

##!@todo this is now the longest-running stage?
##!@todo add option for pre-task slack (in addition to post-task slack).
##! We can't schedule from the back until this is done.
slackEdges = function(schedule, activeWaitConf, critPath, needSlack,
  doubleSlack = F)
{
###!For now, we use the minimum power recorded in any run as active
###!wait power.

###!@todo get a measurement for active wait power and set an active
###wait config for all edge uids

  setkey(schedule, e_uid)
  critEdgeIndices = schedule[J(critPath), which=T]
  ##nonCritEdge_UIDs = setdiff(unique(schedule[, e_uid]), critPath)
  nonCritEdgeIndices = setdiff(1:nrow(schedule), schedule[J(critPath),which=T])

  ## for LP purposes, 0 is equivalent to NA for power. Weight must be
  ## treated differently.
  if(!length(nonCritEdgeIndices) ||
     (!doubleSlack && !missing(needSlack) && length(needSlack) == 0)){
    result = schedule
  } else if(doubleSlack){
    ## doubleSlack forces slack for all tasks
    origEdges = schedule
    slackEdgesPre = data.table::copy(origEdges)
    slackEdgesPost = data.table::copy(origEdges)
    origEdges[, c('src', 'dest', 's_uid', 'd_uid') :=
              list(-e_uid-.5, -e_uid, NA, NA)]
    ##!@todo this will change if we schedule from the back; adjust start times
    slackEdgesPre[, c('e_uid', 'dest', 'd_uid', 'weight') :=
                  list(-e_uid-.5, -e_uid-.5, NA, NA)]
    slackEdgesPost[, c('e_uid', 'src', 's_uid', 'weight', 'start') :=
                   list(-e_uid, -e_uid, NA, NA, start + weight)]
    for(col in names(activeWaitConf)){
      slackEdgesPre[[col]] = activeWaitConf[[col]]
      slackEdgesPost[[col]] = activeWaitConf[[col]]
    }
    result = .rbindlist(list(slackEdgesPre, origEdges, slackEdgesPost))
    rm(slackEdgesPre, slackEdgesPost)
  } else { ## !doubleSlack && length(needSlack) > 0
### !doubleSlack forces slack for non-critical edges in the initial schedule
    origEdges = data.table::copy(schedule[dest %in% needSlack])
    schedule = schedule[!dest %in% needSlack]
    slackEdgesPost = data.table::copy(origEdges)
    origEdges[, c('dest', 'd_uid') := list(-e_uid, NA)]
    slackEdgesPost[, c('e_uid', 'src', 's_uid', 'weight',
                       'start') :=
                   list(-e_uid, -e_uid, NA, NA, start + weight)]
    ## message edges already have conf and power set
    for(col in names(activeWaitConf))
      slackEdgesPost[type=='comp'][[col]] = activeWaitConf[[col]]
    result =
      .rbindlist(list(schedule, origEdges, slackEdgesPost))
  }

  attr(result, 'doubleSlack') <- doubleSlack
  return(result)
}

##!@todo there has to be a better way to do this. It could be
##!parallelized by unique "by" groups.
reduceNoEffect = function(x, x_inv, measurementCols, by, invKey){
  setkeyv(x, invKey)
  setkeyv(x_inv, invKey)
  if('flags' %in% names(x_inv)){
    #xOMP = x[x_inv[bitwAnd(flags, flagBits$omp), list(s_uid)]]
    xOMP = x[OMP_NUM_THREADS > 1]
    xNoOMP = x[OMP_NUM_THREADS == 1]
    #xNoOMP = x[x_inv[!bitwAnd(flags, flagBits$omp), list(s_uid)]]
  } else if('type' %in% names(x)){
    xNoOMP = x[type == 'message']
    xOMP = x[type != 'message']
  }
  #xNoOMP[, OMP_NUM_THREADS:=1]
  cores = getOption('mc.cores')
  if(!is.null(cores) && cores > 1){
    setkeyv(xNoOMP, by)
    u  = unique(xNoOMP[, by, with=F])
    chunks = chunk(u, nrow(u)/cores)
    rm(u)
    xNoOMP =
      .rbindlist(mclapply(chunks, function(ch)
                         ##!@todo fix warning
                         xNoOMP[ch,lapply(.SD, mean), by=by]))
    rm(chunks)
  } else
    xNoOMP = xNoOMP[,lapply(.SD, mean), by=by]
  rbind(xNoOMP, xOMP, use.names=T)
}

..pareto = function(result, convex=T){
  result = result[order(weight, power)][!duplicated(cummin(power))]

  ## this returns a convex power/time frontier, not a convex
  ## power/performance frontier. The first does not imply the second.
  if(nrow(result) > 1 && convex){
    maxWeight = tail(result, 1)[, weight]
    maxPower = max(result[, power])
    result = rbind(result[1], result)
    result[1, c('weight','power') := list(1.01 * maxWeight, 1.01 * maxPower)]

    ## returns convex hull of result with fake point at tail
    result = result[m_chull(weight, power)]
    ## if(debug){
    ##   maxPower = max(result[, power])
    ##   if(tail(result, 1)[, power] != maxPower)
    ##     stop('incorrect assumption\n')
    ## }
    result = head(result, -1)
    if(nrow(result) > 2){
      result[, slope := c(Inf, diff(power)/diff(weight))]
      ## remove neighboring similar slopes
      repeat{
        diffSlope = diff(result[, slope])
        sel = abs(diffSlope) <= 1e-7
        if(any(sel, na.rm=T)){
          result = result[!sel]
        } else
          break
      }
      result[, slope := NULL]
    }
  }
  result
}

pareto = function(edges){
  startTime = Sys.time()

  setkey(edges, e_uid)
  e_uid_count = length(unique(edges[, e_uid]))

  ## return the rows with configurations on the pareto frontier for each edge uid
  .pareto = function(e){
    if(e %% 1000 == 0)
      cat('e_uid', e, 'of', e_uid_count, '\n')
    ## get Pareto frontier; note that this is not necessarily piecewise linear
    result = ..pareto(edges[J(e)])

    ## reduce configurations within .5% of each other
    ##result = result[chull(signif(result[,list(weight,power)], 2))][order(weight, power)]

    
    ## if(is.unsorted(result[, weight]))
    ##   stop('Pareto error: weight ', e)
    ## if(is.unsorted(rev(result[, power])))
    ##   stop('Pareto error: power ', e)
    result
  }
  result = .rbindlist(mclapply(unique(edges[, e_uid]), .pareto))
  cat('Pareto time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  result
}

chunk = function(d, n){
  ##!@todo verify cases where n does not evenly divide d
##  if(inherits(d, 'data.table')){
##    lapply(split(1:nrow(d), ceiling(seq_along(d)/n)), function(s) d[s])
##  } else
  split(d, ceiling(seq_along(d)/n))
}

m_chull = function(x, y){
  x <- cbind(x, y)
  if (nrow(x) == 1)
    return(1L)
  res <- .Call(grDevices:::C_chull, x)
  return(res)
}

loadMergedData = function(filename){
  load(filename, envir=.GlobalEnv)
  setkey(reduced$edges, e_uid)
  setkey(reduced$edges_inv, e_uid)
  setkey(reduced$schedule, e_uid)
}

writeGraphFromSchedule =
  function(sched, critPath, key, compile=T,
           v=sched[,list(vertex=union(src,dest))])
{
  e = sched[,list(e_uid, src,dest,label=paste(rank, '/', e_uid, sep=''))]
  if(!'label' %in% names(v))
    v[,label:=as.character(vertex)]
  else
    v[,label:=as.character(label)]
  
  if(any(sched$e_uid < 0)){
    e$style = 'solid'
    e[e_uid < 0, style := 'dotted']
  }
  
  if(!missing(critPath)){
    e$color = 'black'
    setkey(e, e_uid)
    e[J(critPath), color := 'red']
  }
    
  e$e_uid = NULL
  g=graph.data.frame(e, directed=T,
    vertices=v[, intersect(names(v), c('vertex','label')), with=F])
  key = gsub('[/. ]', '_', key)
  filename = paste(key, '.dot', sep='')
  write.graph(g, filename, format='dot')
  if(compile)
    system(paste('dot -Tpdf -O', filename, sep=' '))
}

plotPerfPower = function(edges, name='Dummy Region', toPDF=T){
  main = paste('Power vs. Performace for', name)
  if(toPDF)
    pdf(paste(gsub('[/. ]', '_', main), '.pdf', sep=''))
  else
    dev.new()
  edges$perf = 1/edges$weight * min(edges$weight)
  plot(edges$power,
       edges$perf,
       xlim=c(min(edges$power)-5, max(edges$power)*1.1),
       ylim=c(0,1),
       main=main,
       xlab='Power (w)',
       ylab='Normalized Performance')
  ep = ..pareto(edges)
  lines(ep$power, ep$perf, col='red', lwd=2)
  threadCounts = unique(edges$OMP_NUM_THREADS)
  edges = edges[order(OMP_NUM_THREADS, cpuFreq)]
  for(tc in threadCounts){
    e = edges[OMP_NUM_THREADS == tc]
    lines(e$power, e$perf)
    text(tail(e$power,1)+4, tail(e$perf,1) - .03, paste(tc, 'threads'))
  }
  ##!@todo key
  if(toPDF)
    dev.off()
}
