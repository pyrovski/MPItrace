powerTime = function(x){
  ranks  = unique(x$rank)
  powerCols = grep('_w$',names(x),value=T)
  cols = c('start','duration',powerCols)
  powerSteps =
    lapply(ranks, function(r){
      x = x[rank == r][order(uid), cols, with=F]
      x = x[,list(start, power=pkg_w+dram_w+pp0_w)]
      list(times = x$start, steps = stepfun(x$start, c(x$power,tail(x$power,1))))
    })
  times = sort(unique(unlist(lapply(powerSteps, '[[', 'times'))))
  steps = lapply(powerSteps, '[[', 'steps')
  powers =
    as.data.table(cbind(start=times,
          power=rowSums(as.data.table(lapply(steps, function(s) s(times))))))
  return(powers)
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
timeslice = function(sched, edges, criticalPath, start=0, length=.01, n=1){
  if(n != 1)
    stop("Not implemented yet")

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

## edges should be in topological order
longest.path = function(edges, vertices, graph){
  if(!inherits(edges, 'data.table'))
    stop('expected data.table')
  if(!inherits(graph, 'igraph'))
    stop('expected igraph')

  vertices$longest = -Inf

  ## via: edge uid
  vertices$via = as.numeric(NA)
  setkey(vertices, vertex)
  vertices[J(1), longest:=0] ## Init
  for(row in 1:nrow(edges)){
    r = edges[row]
##    l = vertices[J(row$dest), longest]
    srcLongest = vertices[J(r$src), longest] + r$weight
    curLongest = vertices[J(r$dest), longest]
    if(srcLongest > curLongest)
      vertices[J(r$dest), via:=r$e_uid]
    vertices[J(r$dest), longest:=max(srcLongest,curLongest)]
  }

  f = function(e){
    via = vertices[J(edges[J(e),src]), via]
    if(!is.na(via))
      return(c(f(via), e))
    else
      e
  }
  
  setkey(edges, e_uid)
  f(vertices[name == 'MPI_Finalize', via])
}

minConf = function(confs){
  ## if we can't simultaneously minimize all knobs, choose the config
  ## with the lowest power. This may be arbitrary for modeled power.
  confs[power > 0][which.min(power)]
}

## orders vertices and edges in topological order and adds start time
## to each edge.
getSchedule = function(edges, vertices=edges[,list(vertex=union(src,dest))]){
  setcolorder(edges,
              c('src','dest',setdiff(names(edges), c('src','dest'))))
  g = graph.data.frame(edges)
  gd = lapply(get.data.frame(g, what='both'), as.data.table)
  gd$vertices$name = as.numeric(gd$vertices$name)
  ts_order = topological.sort(g)

  setkey(vertices)
  vertices = vertices[J(gd$vertices[ts_order])]
  rm(gd)

  edges$start = -Inf
  setkey(edges, src)

  ## define a start time for each edge
  edges[J(vertices$vertex[1]), start:=max(0, start)]
  for(vertex in vertices$vertex){
    outEdges = edges[J(vertex)]
    for(row in 1:nrow(outEdges)){
      startTime = outEdges[row, start + weight]
      edges[J(outEdges[row]$dest), start:=max(startTime, start)]
    }
  }
  edges = edges[J(vertices)]
  return(list(edges=edges, vertices=vertices, graph=g))
}

slackEdges = function(edges, critPath){
  ##!@todo record active wait power in shim library, export to
  ##!glog.dat, and use for slack edges. For now, we use the minimum
  ##!power recorded in any run as active wait power
  slackConfig = minConf(edges[,c(confCols,'power'), with=F])
  slackPower = slackConfig[,power]
  slackConfig = slackConfig[,confCols, with=F]

  setkey(edges, e_uid)
  ##  cols = setdiff(names(edges), confCols)
  nonCritEdges = edges[!J(critPath)]
  critEdges = edges[J(critPath)]

  ## for LP purposes, 0 is equivalent to NA for power. Weight must be
  ## treated differently.
  f = function(row){
    slack = data.table::copy(row)
    row[, c('dest', 'd_uid') := list(-e_uid, NA)]
    if(identical(row[,confCols,with=F], slackConfig)){
      slack[, c('e_uid', 'src', 's_uid', 'weight', 'power') :=
            list(-e_uid, -e_uid, NA, NA, slackPower)]
      rbind(row, slack)
    } else
       row
  }
  nonCritEdges = rbindlist(rowApply(nonCritEdges, f))
  rbind(critEdges, nonCritEdges)
}

##!@todo there has to be a better way to do this
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
  xNoOMP_reduced =
    xNoOMP[,lapply(.SD[, measurementCols, with=F], mean), by=by]
  setkeyv(xNoOMP_unreduced, by)
  setkeyv(xNoOMP_reduced, by)
  xNoOMP = xNoOMP_unreduced[xNoOMP_reduced,, mult='first']
  rbind(xNoOMP, xOMP, use.names=T)
}

pareto = function(edges){
  ##!return the rows with configurations on the pareto frontier for each edge uid
  confs = unique(edges[,confCols, with=F], by=confCols)
  setkey(edges, e_uid)
  f = function(uid_edges){
    uid_edges = uid_edges[order(weight, power)]
    count = 1
    frontier = c(1)
    if(nrow(uid_edges) > 1)
      for(row in 2:nrow(uid_edges))
        if(uid_edges[row, power] < uid_edges[frontier[count], power])
          frontier[count <- count + 1] = row
    uid_edges[frontier]
  }
  edges[, f(.SD), by=e_uid]
}
