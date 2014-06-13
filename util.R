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
  powerSteps = mclapply(ranks, function(r){
    edges = edges[rank == r]
    ##! data.table doesn't order correctly?
    ##[order(start)]
    o = order(edges[,start])
    edges = edges[o]
    times=edges[, start]
    powers = f(edges[, power])
    steps = stepfun(x=times, y=powers)
    list(times=times, steps=steps)
  })
  times = sort(unique(unlist(lapply(powerSteps, '[[', 'times'))))
  steps = lapply(powerSteps, '[[', 'steps')
  powers =
    as.data.table(cbind(start=times,
          power=rowSums(as.data.table(lapply(steps, function(s) s(times))))))

  ## remove runs of identical powers
  cs = cumsum(rle(rev(powers[, power]))$lengths)
  sel = rev(tail(cs, 1) + 1 - cs)
  powers = powers[sel]
  
  return(stepfun(powers[,start], f(powers[, power])))
}

rlePower = function(powers){
  b = rle(powers$power)
  newPowers = powers[1:b$lengths[1]]
  row = 1 + b$lengths[1]
  for(rle_index in 2:length(b$lengths)){
    newPowers = rbind(newPowers, powers[row])
    row = row + b$lengths[rle_index]
  }
  rbind(newPowers, as.data.table(list(start=max(powers$start), power=0)))
  ##newPowers
}

## return sets of indices for each timeslice and how much of each task
## is in each timeslice
timeslice = function(sched, edges, criticalPath, start=0, length=.01, n=Inf){
  sched$end = sched[, start+weight]
  maxStart = max(sched[,start])
  if(start > maxStart)
    stop('Start past last task')
  
  sliceTimes = head(seq(start, maxStart, length), n)
  setkey(sched, e_uid)
  setkey(edges, e_uid)

  f = function(sliceTime){
    nextSlice = sliceTime + length
### There's a package for finding interval overlap, but we only do
### this once.
    int =
      sched[start >= sliceTime & end <= nextSlice,
            list(e_uid)] ## interior
    int[,c('left', 'right') := list(0, 1)]
    ext =
      sched[start <  sliceTime & end >  nextSlice,
            list(e_uid, left=(sliceTime-start)/weight,
                 right=(end-nextSlice)/weight)] ## exterior
    ## left overlap
    left =
      sched[start <  sliceTime & end >   sliceTime & end <= nextSlice,
            list(e_uid, left=(sliceTime-start)/weight)]
    left[, right := 1]
    ## right overlap
    right =
      sched[start >= sliceTime & start < nextSlice & end >  nextSlice,
            list(e_uid, right=(nextSlice-start)/weight)]
    right[, left := 0]
    setcolorder(right, names(int))
    
    result = rbind(int, ext, left, right)

###!we need to modify all edges according to the fraction in the
###!timeslice, not just the scheduled edges. So, find out which edges
###!are in the timeslice, then determine their splits. Possible splits
###!depend on the type of overlap.

    return(edges[result][,c('weight', 'left', 'right') :=
                         list((right-left)*weight, NULL, NULL)])
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
  g = graph.data.frame(edges)
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
  vertices$start = -Inf
  vertices[J(1), start := 0]
  vertices[, via:=as.numeric(NA)]
  count = 0
  progress = 0
  
###!@todo this could be faster if we selected only the relevant set of
###!vertices for each src, then merged after the loop
  for(vertex in vertices_TO$vertex){
    ##!@todo some of these columns are not used
    outEdges = edges[J(vertex), list(src, dest, e_uid, s_uid, d_uid, type, weight)]
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
    vertices[J(v[start_u > start,
                 list(start_u, e_uid),
                 keyby=vertex]),
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
    warning('Some edges not assigned a start time!\n', immediate. = T)
  
  if(doCritPath){
    cat('Critical path\n')

###!igraph doesn't support finding the shortest path with negative
###!edge weights.

    setkey(edges, e_uid)
    c_vertex = 2
    critPath = c()
    while(c_vertex != 1){
      via = vertices[J(c_vertex), via]
      c_vertex = edges[J(via), src, mult='first']
      critPath = c(via, critPath)
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
  ##!return the rows with configurations on the pareto frontier for each edge uid
  f = function(uid_edges){
    uid_edges = uid_edges[order(weight, power)]

    frontier = NULL
    while(nrow(uid_edges) > 0){
      frontier = rbind(frontier, uid_edges[1])
      uid_edges = uid_edges[power < tail(frontier, 1)[, power]]
    }
    frontier
  }
  setkey(edges, e_uid)
  rbindlist(mclapply(unique(edges[, e_uid]), function(e) f(edges[J(e)])))
}

chunk = function(d, n){
  ##!@todo verify cases where n does not evenly divide d
  split(d, ceiling(seq_along(d)/n))
}
