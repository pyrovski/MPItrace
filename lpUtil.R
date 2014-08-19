require('igraph')
require('rjson')
require('data.table')
require('parallel')
source('./util.R')
source('~/local/bin/pbutils.R')

readLP = function(filename){
  a = fromJSON(file=filename)

  f = function(b){
    arraySel = grep('[[]', names(b))
    singleSel = setdiff(1:length(b), arraySel)
    if(length(arraySel)){
      arrayNames = gsub('[[].*[]]', '', names(b)[arraySel])
      uArrayNames = unique(arrayNames)
      arrayVars = nnapply(uArrayNames, function(name){
        b = b[grep(paste(name, '[[]', sep=''), names(b))]
        map = strsplit(gsub('[]]', '', names(b)), '[[]')
        indices = as.numeric(sapply(map, '[[', 2))
        b = rbindlist(lapply(b, as.data.table))
        b$index = indices
        if('Id' %in% names(b))
          b[, Id := NULL]
        b
      })
    } else
      arrayVars = NULL
    if(length(singleSel)){
      singleVars = lapply(b[singleSel], function(e){
        b = as.data.table(e)
        if('Id' %in% names(b))
          b[, Id := NULL]
        if(ncol(b) == 1 && identical(names(b), 'Value'))
          b = b[[1]]
        b
      })
    } else
      singleVars = NULL
    c(arrayVars, singleVars)
  }
  a$Solution[[2]]$Variable = f(a$Solution[[2]]$Variable)
  a$Solution[[2]]$Constraint = f(a$Solution[[2]]$Constraint)
  a$Solution[[2]]$Objective = f(a$Solution[[2]]$Objective)
  a
}

reconcileLP = function(resultFile, timesliceFile, keepAll=F){
  result = readLP(resultFile)
  vertexStartTimes = result$Solution[[2]]$Variable$vertexStartTime
  result$Solution[[2]]$Variable$vertexStartTime = NULL
  setnames(vertexStartTimes, c('index', 'Value'), c('vertex', 'start'))
  setkey(vertexStartTimes, vertex)

  tryCatch(load(timesliceFile), finally=NULL)
  vertices = slice[, list(vertex=union(src, dest))]
  setkey(vertices, vertex)
  vertices = vertexStartTimes[vertices]
  rm(vertexStartTimes)
  ##!@todo this is incorrect; some vertices are not present in this
  ##!timeslice or the next; just because they're destinations doesn't
  ##!mean we need to assign them start times.
  vertices[is.na(start), start := 0.0]
  
  unconstrained = vertices[start > .95]
  if(nrow(unconstrained) > 0){
      cat('unconstrained vertices:', resultFile, '!\n')
      print(unconstrained)
  }
  rm(unconstrained)

  e_uids = unique(slice[, e_uid])
  
  taskDuration = result$Solution[[2]]$Variable$taskDuration
  setnames(taskDuration, c('index', 'Value'), c('e_uid', 'lpWeight'))
  setkey(taskDuration, e_uid)
  taskDuration = taskDuration[J(e_uids)]
  taskDuration[is.na(lpWeight), lpWeight := 0]

  taskPower = result$Solution[[2]]$Variable$taskPower
  setnames(taskPower, c('index', 'Value'), c('e_uid', 'lpPower'))
  setkey(taskPower, e_uid)
  taskPower = taskPower[J(e_uids)]
  setkey(taskPower, e_uid)
  ## these should not exist, at least for comp edges
  ##!@todo warn on NA power for comp edges
  taskPower = taskPower[slice[,head(.SD, 1),keyby=e_uid,.SDcols=c('type')]]
  if(nrow(taskPower[is.na(lpPower) & type == 'comp']) > 0){
    stop('LP should provide all comp task power entries')
  }
  taskPower[is.na(lpPower) & type == 'comp', lpPower := 0]
  taskPower[, type := NULL]

  setkey(taskDuration, e_uid)
  task = merge(taskDuration, taskPower, all=T)
  rm(taskPower, taskDuration)
  
  setkey(slice, e_uid, weight, power)
  setkey(task, e_uid, lpWeight, lpPower)
  f = function(a, b) abs(a-b) < 1e-8
  edges = lapply(e_uids, function(u){
    s = slice[J(u)]
    if(nrow(s) == 1){
      s[, frac := 1]
      return(s)
    }
    lp = task[J(u)]

    if(keepAll){
      return(list(slice=slice, lp=lp))
    }
    
    unconstrained = lp[lpWeight > .9 & (lpWeight %% 1) < .1]
    if(nrow(unconstrained) > 0){
      cat('unconstrained weight(s)!\n')
      print(unconstrained)
    }
    
    ##!@todo this can be done with multiple e_uids at once

    ##! this needs to be approximate
    ##m = s[lp, nomatch=0]
    m = s[f(weight, lp$lpWeight) & f(power, lp$lpPower)]
    if(nrow(m) > 0){
      m[, frac := 1]
      return(m)
    }

    ##!@todo figure out how to get Pyomo to be more precise with its output

    ##! can re-adjust lp weight based on selected power
    m = rbind(head(s[power < lp$lpPower], 1), tail(s[power > lp$lpPower], 1))
    fastFrac = (lp$lpPower - m[1, power])/diff(m[, power])
    slowFrac = 1 - fastFrac
    m$frac = c(fastFrac, slowFrac)
    ##! adjust weight by frac
    m[, weight := weight * frac]
###! m should contain two rows; one for each configuration neighboring
###! the LP-selected power/performance point
    return(m)
  })
  if(!keepAll)
    edges = rbindlist(edges)
  
  list(
    ##result=result,
    ##slice=slice,
    vertices=vertices[order(start)],
    edges=edges)
}


timeStr = '[0-9]+[.][0-9]+'

##!@todo save results from this function, check for newer inputs than previous result
readCommandResults = function(command, ...){
  cat(command, '\n')
  resultFiles =
    list.files(pattern=
               paste(command, '_', timeStr, '[.]p.*w[.]results$', sep=''))
  powerLimits =
    unique(sub('w[.]results$', '',
               sub(paste(command, '_', timeStr, '[.]p', sep=''), '', resultFiles)))
  prefixes = unique(sub('[.]p.*w[.]results$', '', resultFiles))
  timesliceFiles = paste(prefixes, '.Rsave', sep='')
  times = sub(paste(command, '_', sep=''), '', prefixes)
  f = function(powerLimit){
    cat(powerLimit, 'w', '\n')
    resultFiles =
      list.files(pattern=
                 paste(command, '_', timeStr, '[.]p', powerLimit, 'w[.]results$',
                       sep=''))
    times = sub('[.]p.*w[.]results$', '',
      sub(paste(command, '_', sep=''), '', resultFiles))
    result = mcmapply(reconcileLP, resultFiles, timesliceFiles, ..., SIMPLIFY=F)
    names(result) = times
    result
  }
  nnapply(powerLimits, f)
}

lpGo = function(...){
  load('../mergedEntries.Rsave', envir=.GlobalEnv)
  files = list.files(pattern='.*[.]results$')
  commands <<- unique(sub(paste('_', timeStr, '[.]p.*w[.]results', sep=''),'',files))
  results <<- nnapply(commands, readCommandResults, ...)
  NULL
}

lpMerge = function(slices, name){
  edges =
    rbindlist(napply(slices, function(e, name) {
      e$edges$ts = name
      e$edges
    }, mc=T))
  vertices =
    rbindlist(napply(slices, function(e, name) {
      e$vertices$ts =name
      e$vertices
    }, mc=T))
  tsDuration = vertices[, .SD[which.max(start)], by=ts]
  tsDuration[, vertex := NULL]
  setnames(tsDuration, 'start', 'tsEnd')
  setkey(tsDuration, ts)
  tsDuration[, tsEnd := cumsum(tsEnd)]
  tsDuration$tsStart = 0
  tsDuration$tsStart[2:nrow(tsDuration)] = head(tsDuration[, tsEnd], -1)

  setkey(vertices, ts)
  vertices = vertices[tsDuration[, list(ts, tsStart)]]
  vertices[, c('start', 'tsStart') := list(start + tsStart, NULL)]

  setnames(vertices, 'vertex', 'src')
  setkey(vertices, ts, src)
  setkey(edges, ts, src)
  edges = vertices[edges]
  setnames(vertices, 'src', 'vertex')

  ##!@todo renumber vertices across timeslices
  edges[, c('src', 'dest') := list(as.character(src), as.character(dest))]
  edges[, ts := as.character(.GRP), by=ts]
  edges[, ts := as.integer(ts)]
  edges[splitDest == T, dest := paste(src, '_', ts, 's', sep='')]
  edges[splitSrc == T, src := paste(src, '_', ts-1, 's', sep='')]
  edges[, c('splitSrc', 'splitDest') := list(NULL, NULL)]

  ## just to be consistent
  ##edges[splitSrc==F, src := paste(src, ts, sep='_')]
  ##edges[splitDest==F, dest := paste(dest, ts, sep='_')]
  
  ## assign new vertices to split-config edges from each timeslice,
  ## rename edge uids to be unique across timeslices
  edges = edges[order(ts, e_uid, -frac)]
  edges[, second := F]
  edges =
    edges[,if(.N ==2){
      e = copy(.SD)
      e[1, c('dest') := list(paste(src, '.', sep=''))]
      e[2, c('src', 'start', 'second') :=
        list(paste(src, '.', sep=''), start + e[1, weight], T)]
      e
    } else {
      .SD
    }, by=list(e_uid, ts)]
  edges[, e_uid := as.character(e_uid)]
  edges[, e_uid := paste(ts, e_uid, sep='_')]
  edges[second == T, e_uid := paste(e_uid, '.', sep='')]
  
  ##!@todo assign weights to slack edges

  vertices = edges[, list(vertex=union(src, dest))]
  setkey(vertices, vertex)
  setkey(edges, src)
  vertices[edges[,list(src, start)], start:=start]
  vertices[J('2'), start:=edges[dest=='2', max(start+weight)]]

  pt = powerTime(edges, vertices)
  plotPowerTime(pt, name=name)
  
  return(list(edges = edges,
              vertices = vertices,
              pt = pt))
}

if(!interactive()){
  lpGo()
  results = lapply(results, function(e) lapply(e, lpMerge))
}
