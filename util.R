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
    result = rbindlist(mclapply(sets, function(s) x[s]))
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

powerTime = function(edges){
  ranks = unique(edges[, rank])
  f = function(x) c(x, 0)
  rf = function(r){
    edges = edges[rank == r & power > 0]
    ##! data.table doesn't order correctly?
    ##[order(start)]
    o = order(edges[,start])
    edges = edges[o]
    times=edges[, start]
    powers = f(edges[, power])
    steps = stepfun(x=times, y=powers)
    list(times=times,
         ##powers=powers,
         steps=steps)
  }
  powerSteps = mclapply(ranks, rf)
  times = sort(unique(unlist(lapply(powerSteps, '[[', 'times'))))
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

powerStats = function(edges){
### need to assign start times for each set

  ##fastest
  sched = getSchedule(edges[,.SD[which.min(weight)],by=e_uid])
  sched = slackEdges(sched$edges, sched$critPath)
  powersMinTime = powerTime(sched)
  
  ## max power
  sched = getSchedule(edges[,.SD[which.max(power)],by=e_uid])
  sched = slackEdges(sched$edges, sched$critPath)
  powersMaxPower = powerTime(sched)

  ## min power
  sched = getSchedule(edges[,.SD[which.min(power)],by=e_uid])
  sched = slackEdges(sched$edges, sched$critPath)
  powersMinPower = powerTime(sched)

  list(minTime  = powersMinTime,
       maxPower = powersMaxPower,
       minPower = powersMinPower)
}

## return sets of indices for each timeslice and how much of each task
## is in each timeslice. Slack edges have zero weight, so they will
## never be bisected by a timeslice. They are important to include
## because they contribute to power consumption.

##!@todo for the following situation, make sure to include the
##!preceding slack edge in the second timeslice:
##|        x|xx      x|x          
##|   xxxxx |  xxxxxx |            
timeslice = function(sched, vertices, edges, criticalPath,
  ts_start=0, length=.01, n=Inf){
  setkey(vertices, vertex)
  setkey(sched, dest)
  sched = sched[vertices[, list(vertex, deadline=start)]]
  if(nrow(sched[start > deadline]) > 0){
    stop('invalid deadline!\n')
  }
  sched[, wslack:=deadline-start]
  maxStart = max(sched[,start])
  if(ts_start > maxStart)
    stop('Start past last task')
  
  sliceTimes = head(seq(ts_start, maxStart, length), n)
  setkey(sched, e_uid)
  setkey(edges, e_uid)

  f = function(sliceTime){
    nextSlice = sliceTime + length
### There's a package for finding interval overlap, but we only do
### this once.

    ## interior edges: entirely included in the interval
    int =
      sched[start >= sliceTime & deadline <= nextSlice,
            list(e_uid, weight, wslack, right=wslack)]
    int[, left := 0]

    ## exterior edges: edges include interval
    ext =
      sched[start <  sliceTime & deadline >  nextSlice,
            list(e_uid, weight, wslack, left=sliceTime-start,
                 right=nextSlice-start)]

    ## left overlap
    left =
      sched[start <  sliceTime & deadline >   sliceTime & deadline <= nextSlice,
            list(e_uid, weight, wslack, left=sliceTime-start, right=wslack)]
    
    ## right overlap
    right =
      sched[start >= sliceTime & start < nextSlice & deadline >  nextSlice,
            list(e_uid, weight, wslack, right=nextSlice-start)]
    right[, left := 0]
    
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

    result[weight != 0, c('frac','left','right') :=
           list((right-left)/wslack, NULL, NULL)]
    result[weight == 0, frac := 0]
    result[, c('weight', 'wslack') := list(NULL, NULL)]
    result = edges[result][,c('weight','frac') := list(frac*weight, NULL)]
    if(nrow(result[weight < 0]) > 0)
      stop('Negative-weight edge(s)\n')
    return(result)
  }

  result = nnapply(sliceTimes, f)
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
  
  cat('Graph construction\n')
  edges = data.table::copy(edges)
  setcolorder(edges,
              c('src','dest',setdiff(names(edges), c('src','dest'))))
  if(any(edges[, weight] < 0)){
    stop('Negative-weight edge(s)!')
  }
  g =
    graph.data.frame(edges[,c('src','dest','weight','power','e_uid',
                              's_uid','d_uid','rank','type',
                              confCols), with=F])
  gd = lapply(get.data.frame(g, what='both'), as.data.table)
  gd$vertices$name = as.numeric(gd$vertices$name)
  cat('Topological ordering\n')
  ts_order = topological.sort(g)
  if(no.clusters(g) > 1)
    stop('Graph has more than one cluster!', immediate. = T)
  rm(g)
  
  setkey(vertices)
  vertices_TO = data.table::copy(vertices[J(gd$vertices[ts_order])])
  rm(gd)

  setkey(edges, src)

  ## define a start time for each edge
  cat('Start times\n')
  startTime = Sys.time()

  vertices$start = -Inf
  vertices[J(1), start := 0]
  vertices[, via:=as.integer(NA)]
  count = 0
  progress = 0
  
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
    ##!@todo some of these columns are not used
    outEdges = edges[J(vertex), list(src, dest, e_uid, weight)]
    if(nrow(outEdges) < 1)
      next
    setkey(outEdges, src)

    ## get start times for src vertices
    v = vertices[outEdges]

    ## get potential start times for dest vertices
    v[, start_u:=start+weight]

    ## resolve cases with multiple edges between src and dest
    v = v[, .SD[which.max(start_u)], keyby=dest]
    ## should now have one edge per dest in v (already have a single src)

    ## update start times for dest vertices
    v = vertices[v[, list(dest, start_u, e_uid)]]
    ## join on the dest vertex from v and 'vertex' from vertices
    
    vertices[J(v[start_u > start, list(start_u, e_uid), keyby=vertex]),
             c('start', 'via') := list(start_u, e_uid)]
    
    if(!count %% 1000){
      n_progress = round(count/nrow(vertices)*100)
      if(progress != n_progress)
        cat(n_progress, '%\n')
      progress = n_progress
    }
    count = count + 1
  }
  cat(100, '%\n')

###!@todo this assumes edges start as soon as their dependencies are
###!satisfied. If we want to schedule from the back, edges should
###!start as late as possible; for each edge:
###! start = vertices[J(dest), start] - weight
  edges = edges[J(vertices[, list(vertex, start)])]
 
  if(any(edges$start == -Inf))
    stop('Some edges not assigned a start time!\n')
  cat('Start times time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')

  if(doCritPath){
    cat('Critical path\n')

###!igraph doesn't support finding the shortest path with negative
###!edge weights.

    setkey(edges, e_uid)
    c_vertex = 2
    critPath = c()
    while(c_vertex != 1){
      v = vertices[J(c_vertex)][, via]
      c_vertex = edges[J(v), src, mult='first']
      critPath = c(v, critPath)
    }
  } else
    critPath=NULL

  return(list(edges=edges, vertices=vertices, critPath=critPath))
}

slackEdges = function(edges, critPath){
###!For now, we use the minimum power recorded in any run as active
###!wait power

  setkey(edges, e_uid)
  ##  cols = setdiff(names(edges), confCols)
  ##nonCritEdges = edges[!J(critPath)]
  nonCritEdge_UIDs = setdiff(unique(edges[, e_uid]), critPath)
  critEdges = edges[J(critPath)]

  ## for LP purposes, 0 is equivalent to NA for power. Weight must be
  ## treated differently.
  f = function(u){
    confs = edges[J(u)]
    slack = confs[which.min(power)]
    confs[, c('dest', 'd_uid') := list(-e_uid, NA)]
    slack[, c('e_uid', 'src', 's_uid', 'weight', 'start') :=
          list(-e_uid, -e_uid, NA, NA, start + weight)]
    rbind(confs, slack)
  }
  if(length(nonCritEdge_UIDs)){
    nonCritEdges = rbindlist(mclapply(nonCritEdge_UIDs, f))
    result = rbind(critEdges, nonCritEdges)
  } else
    result = critEdges
  return(result)
}

##!@todo there has to be a better way to do this. It could be
##!parallelized by unique "by" groups.
reduceNoEffect = function(x, measurementCols, nonMeasurementCols, by){
  if('flags' %in% names(x)){
    xOMP = x[flags & flagBits$omp]
    xNoOMP = x[!flags & flagBits$omp]
  } else if('type' %in% names(x)){
    xNoOMP = x[type == 'message']
    xOMP = x[type != 'message']
  }
  xNoOMP[, OMP_NUM_THREADS:=1]
  xNoOMP_unreduced = xNoOMP[, nonMeasurementCols, with=F]
  xNoOMP = xNoOMP[, c(measurementCols, by), with=F]
  cores = getOption('mc.cores')
  if(!is.null(cores) && cores > 1){
    setkeyv(xNoOMP, by)
    u  = unique(xNoOMP[, by, with=F])
    chunks = chunk(u, nrow(u)/cores)
    rm(u)
    xNoOMP_reduced =
      rbindlist(mclapply(chunks, function(ch)
                         xNoOMP[ch,lapply(.SD, mean), by=by]))
    rm(chunks)
  } else
    xNoOMP_reduced = xNoOMP[,lapply(.SD, mean), by=by]
  setkeyv(xNoOMP_unreduced, by)
  setkeyv(xNoOMP_reduced, by)
  xNoOMP = xNoOMP_unreduced[xNoOMP_reduced, mult='first']
  rbind(xNoOMP, xOMP, use.names=T)
}


pareto = function(edges){
  startTime = Sys.time()

  setkey(edges, e_uid)
  e_uid_count = length(unique(edges[, e_uid]))

  ##!return the rows with configurations on the pareto frontier for each edge uid
  .pareto = function(e){
    if(e %% 1000 == 0)
      cat('e_uid', e, 'of', e_uid_count, '\n')
    result = edges[J(e)][order(weight, power)][!duplicated(cummin(power))]

    ##!@todo reduce configurations within .5% of each other
    ##result = result[chull(signif(result[,list(weight,power)], 2))][order(weight, power)]

    if(nrow(result) > 1){
      ##slopes = diff(result[, power])/diff(result[, weight])
      result = result[order(weight, power)][!duplicated(cummax(diff(power)/diff(weight)))]
    }

    ## if(is.unsorted(result[, weight]))
    ##   stop('Pareto error: weight ', e)
    ## if(is.unsorted(rev(result[, power])))
    ##   stop('Pareto error: power ', e)
    result
  }
  result = rbindlist(mclapply(unique(edges[, e_uid]), .pareto))
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
