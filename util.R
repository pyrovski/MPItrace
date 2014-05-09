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
timeslice = function(sched, start=0, length=.01, n=1){
  if(n != 1)
    stop("Not implemented yet")

  sched$end = sched[, start+weight]
  maxStart = max(sched[,start])
  if(start > maxStart)
    stop('Start past last task')
  
  sliceTimes = seq(start, maxStart, length)

  f = function(sliceTime){
    nextSlice = sliceTime + length
### There's a package for finding interval overlap, but we only do this once.
    int = sched[start >= sliceTime & end <= nextSlice, which=T] ## interior
    ext = sched[start <  sliceTime & end >  nextSlice, which=T] ## exterior
    ## left overlap
    left = sched[start <  sliceTime & end >   sliceTime & end <= nextSlice,
      which=T]
    ## right overlap
    right = sched[start >= sliceTime & start < nextSlice & end >  nextSlice,
      which=T]

    result = list(int=int, ext=ext, left=left, right=right)
    
###!@todo define overlap fractions (depends on overlap type). For
###!dependent timeslices, we don't care about the right side, just the
###!left side. We only get to decide configurations for tasks that
###!start in the current timeslice (interior and right types); other
###!tasks retain their previous configurations.
    result
  }
  
  result = lapply(sliceTimes, f)
  result
}
