#!/usr/bin/env Rscript

require('igraph')
require('rjson')
require('data.table')
require('parallel')
source('./util.R')
source('./read.R')
source('~/local/bin/pbutils.R')

ilpFileTypes = c('duration', 'edges')
fixedLPFileTypes = c('duration', 'edges')

.solveLP = function(entrySpaceRow, nodePowerLimits, forceLocal=F, fixedOnly=T, roundMode='step'){
  cuts =
    sub('.edges.csv','',
        Sys.glob(paste(gsub('[/. ]', '_', entrySpaceRow$key),
                       '*.edges.csv', sep='')))
  ilpCuts = grep('ILP', cuts, v=T)
  fixedCuts = setdiff(cuts, ilpCuts)
  rm(cuts)
  mclapply(fixedCuts, function(cut){
    mclapply(nodePowerLimits, function(pl) {
      edgesFile = paste(cut, '.p', format(pl, nsmall=1), 'w.edges', sep='')
      cutEdges = paste(cut, '.edges.csv', sep='')
      if(file.exists(edgesFile) && file.info(edgesFile)$mtime > file.info(cutEdges)$mtime){
        cat(edgesFile, 'exists\n')
        return(NULL)
      }
      command =
        paste('prefix=', cut, ' powerLimit=', pl*entrySpaceRow$ranks,
              if(forceLocal) ' FORCELOCAL=1',
              ' roundMode=', roundMode,
              ' ./fixed.sh', sep='')
      print(command)
      system(command, intern=T)
    })
  })
}

solveLP = function(...){
  load('../mergedEntries.Rsave')
  mcrowApply(entrySpace, .solveLP, ...)
}

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
        b = .rbindlist(lapply(b, as.data.table))
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

#!@todo adapt for flow ilp, fixed ilp
reconcileLP = function(resultFile, timesliceFile, powerLimit, mode='split'){
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

  result = list()
  
  if(mode=='keepAll'){
    setkey(slice, e_uid)
    setkey(task, e_uid)
    ##edges = slice[task]
    result$edges = slice
    result$lp = task
  } else { ## mode != 'keepAll'
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
      if(mode=='combined'){
        ## find a single config that is closest to the LP
        m[, dist := sqrt(((power-lp$lpPower)/lp$lpPower)^2+((weight - lp$lpWeight)/lp$lpWeight)^2)]
        m = m[which.min(dist)]
        m$dist = NULL
        m$frac=1
        return(m)
      } else if(mode == 'combinedLE'){
        ## find a single config that is always under the power constraint
        m = m[power <= powerLimit, .SD[which.min(weight)], by=e_uid]
        m$frac=1
        return(m)
      } else if(mode == 'split'){  ## split configs
        fastFrac = (lp$lpPower - m[1, power])/diff(m[, power])
        slowFrac = 1 - fastFrac
        m$frac = c(fastFrac, slowFrac)
        ##! adjust weight by frac
        m[, weight := weight * frac]
###! m should contain two rows; one for each configuration neighboring
###! the LP-selected power/performance point
        return(m)
      }
    })
    edges = .rbindlist(edges)
    result$edges = edges
  }

  result$vertices=vertices[order(start)]
  result
}


timeStr = '[0-9]+[.][0-9]+'

##!@todo save results from this function, check for newer inputs than previous result
##!@todo make sure result files are newer than csv and Rsave inputs
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
    powerLimit = as.numeric(powerLimit)
    cat(powerLimit, 'w', '\n')
    resultFiles =
      list.files(pattern=
                 paste(command, '_', timeStr, '[.]p', powerLimit, 'w[.]results$',
                       sep=''))
    times = sub('[.]p.*w[.]results$', '',
      sub(paste(command, '_', sep=''), '', resultFiles))
    result = mcmapply(reconcileLP, resultFiles, timesliceFiles, powerLimit, ..., SIMPLIFY=F)
    names(result) = times
    result
  }
  nnapply(powerLimits, f)
}

lpGo = function(...){
  load('../mergedEntries.Rsave', envir=.GlobalEnv)
  files = list.files(pattern='.*[.]results$')
  commands <<- unique(sub(paste('_', timeStr, '[.]p.*w[.]results', sep=''),'',files))
  nnapply(commands, readCommandResults, ...)
}

ilpGo = function(pattern='.*', powerLimitMin=0, ...){
  ## get all files, then filter by prefix, then by power limit and cut
  load('../mergedEntries.Rsave', envir=.GlobalEnv)
  cutPattern = 'cut_[0-9]+'
  plPattern = 'p.*w'
  pattern = paste(pattern, '.*[.]duration$', sep='')
  files = list.files(pattern=pattern)
  ##duration'
  prefixes =
    unique(sub(paste('(.*)', cutPattern, plPattern, 'duration', sep='[.]'),
               '\\1', files))

  if(!length(prefixes)){
    warning('no ILP result files!', immediate.=T)
    return(NULL)
  }
  
  nnapply(prefixes, function(prefix){
    fixed = length(grep('fixedLP', prefix)) > 0
    
    files = list.files(pattern=paste(prefix, cutPattern, plPattern,
                         'duration$', sep='[.]'))
    powerLimits =
      sub('p([0-9.]+)w', '\\1',
          unique(sub(paste('.*', cutPattern, paste('(', plPattern, ')', sep=''),
                           'duration$', sep='[.]'),
                     '\\1', files)))

    ##!@todo this needs to be modified to reformat the power limits with a trailing zero
    if(powerLimitMin > 0){
      powerLimitFloat = as.numeric(powerLimits)
      cat('ignoring power limits: ', powerLimits[powerLimitFloat < powerLimitMin], '\n')
      powerLimits = powerLimits[powerLimitFloat >= powerLimitMin]
    }
    
    cuts =
      sort(as.numeric(
        sub('cut_([0-9]+)', '\\1',
            unique(sub(paste(prefix, paste('(',cutPattern, ')', sep=''),
                             plPattern, 'duration$', sep='[.]'),
                       '\\1', files)))))

    expectedCuts = as.integer(read.table(paste(prefix, '.cuts.csv', sep=''), h=F)[[1]])
    
    nnapply(powerLimits, function(powerLimit){
      plPattern = paste('p', powerLimit, 'w', sep='')
      presentCuts =
        list.files(pattern=paste(prefix, cutPattern,
                     plPattern,
                     'edges$', sep='[.]'))
      presentCuts = as.integer(sub('.*cut_([0-9]+).*', '\\1', presentCuts))
      
      if(length(setdiff(expectedCuts, presentCuts))){
        errMsg =
          paste(prefix, '@', powerLimit, 'w:\nmissing cuts!\n',
                paste(setdiff(expectedCuts, presentCuts), collapse=' '), sep='')
        stop(errMsg)
      }
      
      nnapply(cuts,
              function(cut){
                if(fixed)
                  fileTypes = fixedLPFileTypes
                else
                  fileTypes = ilpFileTypes
                
                nnapply(fileTypes, function(fileType){
                  filename =
                    paste(prefix, paste('cut_', cut, sep=''),
                          paste('p', powerLimit, 'w', sep=''),
                          fileType, sep='.')
                  tryCatch(
                    as.data.table(
                      read.table(filename, h=T, sep=',', strip.white=T)),
                    error=function(e){
                      warning('failed to read ', filename, immediate.=T)
                      NULL
                    }, finally=NULL)
                }
                        )
              },
              mc=T)
    }
            )
  }
          )
}

##!@todo this function assumes that we don't alter the schedule from
##!the LP.  For modes other than the default, this may not be true,
##!and we need to recompute start times and slack edges.
lpMerge = function(slices, name){
  edges =
    .rbindlist(napply(slices, function(e, name) {
      e$edges$ts = name
      e$edges
    }, mc=T))
  
  vertices =
    .rbindlist(napply(slices, function(e, name) {
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

  ## renumber vertices across timeslices
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
  if('frac' %in% names(edges))
    edges = edges[order(ts, e_uid, -frac)]
  else
    edges = edges[order(ts, e_uid)]
  
  edges[, second := F]
  edges[, orig_e_uid := e_uid]
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
  e_uid_map = data.table(orig=edges$e_uid)
  edges[, e_uid := as.character(e_uid)]
  edges[, e_uid := paste(ts, e_uid, sep='_')]
  edges[second == T, e_uid := paste(e_uid, '.', sep='')]
  e_uid_map$new = edges$e_uid
  
  ##!@todo assign weights to slack edges

  vertices = edges[, list(vertex=union(src, dest))]
  setkey(vertices, vertex)
  setkey(edges, src)
  vertices = vertices[unique(edges[,list(src, start)])]
  vertices =
    rbind(vertices, data.table(vertex='2',
                               start=edges[dest=='2', max(start+weight)]))

  pt = powerTime(edges, vertices)
  plotPowerTime(pt, name=name)
  
  return(list(edges = edges,
              vertices = vertices,
              pt = pt))
}

lpMergeAll = function(commands_powers){
  napply(commands_powers,
         function(x, name){
           command = name
           napply(x,
                  function(x, name){
                    powerLimit = name
                    lpMerge(x, name=paste(command, powerLimit))
                  }, mc=T)
         })
}

loadAndMergeLP = function(){
  results <<- lpGo()
  resultsMerged <<- lpMergeAll(results)
  resultsOneConf <<- lpGo(mode='combined')
  resultsMergedOneConf <<- lpMergeAll(resultsOneConf)
  resultsOneConfLE <<- lpGo(mode='combinedLE')
  resultsMergedOneConfLE <<- lpMergeAll(resultsOneConfLE)
}

# note: this also merges fixedLP results. I'm lazy.
loadAndMergeILP = function(...){
  resultsILP = ilpGo(...)

  ## retain only complete cuts
  ##!@todo get list of expected cuts, warn if any missing
  
  
  f = function(x, depth){
    if(depth == 1){
      if(any(sapply(x, is.null)))
        NULL
      else
        x
    }
    else {
      result = lapply(x, f, depth-1)
      result = result[!sapply(result, is.null)]
      if(length(result))
        result
      else
        NULL
    }
  }
  resultsILP = f(resultsILP, 4)

  fixed = grep('fixedLP', names(resultsILP))
  ilp = setdiff(seq(length(resultsILP)), fixed)

  resultsFixedLP = resultsILP[fixed]
  resultsILP = resultsILP[ilp]

  f = function(cuts, fileTypes) 
    nnapply(
      fileTypes,
      function(fileType)
      .rbindlist(
        napply(
          cuts,
          function(x, name){
            result=x[[fileType]]
            if(is.null(result)) return(result)
            result[, cut:=as.numeric(name)]
            result
          }
          )
        )
      )
  
  ## merge cuts
  resultsILPMerged = lapply(resultsILP, lapply, f, ilpFileTypes)
  resultsFixedLPMerged = lapply(resultsFixedLP, lapply, f, fixedLPFileTypes)

  ## propagate event times across cuts.

  ##!@todo collectives are numbered in order of their occurrence, but
  ##!e.g. MPI_Waitall()s may not be. As the cut names are vertex
  ##!labels, we can get the vertex ordering from mergedData.
##   f = function(x){
## ### if we somehow avoid zero-length slack edges, the end time of a cut
## ### will not correspond to the start time of its last vertex
##     ## place event start times within each cut

##     setkey(x$duration, cut)
##     x$duration[, cutEnd:=cumsum(duration)]
##     x$duration[, cutStart:=c(0, head(cutEnd, -1))]
##     setkey(x$events, cut, event)
##     x$events = x$duration[x$events, list(event, cut, start=start+cutStart)]
##     setkey(x$edges, cut, event)
##     setkey(x$events, cut, event)
##     x$edges = x$events[x$edges]
##     ## renumber events, remove cut column
##     x$activeEvents = x$edges[, list(event=unique(event)), by=cut]
##     setkey(x$activeEvents, cut, event)
##     x$activeEvents = x$activeEvents[, list(newEvent=.GRP-1), by=list(cut, event)]
##     x$edges[x$activeEvents, c('event', 'cut') := list(newEvent, NULL)]
##     x$events = x$events[x$activeEvents]
##     x
##   }
## #  resultsILPMerged <<- lapply(resultsILPMerged, lapply, f)

  assign('resultsILPMerged', resultsILPMerged, envir=.GlobalEnv)
  assign('resultsFixedLPMerged', resultsFixedLPMerged, envir=.GlobalEnv)
  NULL
}

accumulateCutStarts = function(x, orderedCuts){
  setkey(x$duration, cut)
  x$duration = x$duration[J(orderedCuts)]
  x$duration$cutEnd = as.numeric(NA)
  x$duration$cutStart = as.numeric(NA)
  x$duration[duration != 'infeasible', cutEnd:=cumsum(duration)]
  x$duration[duration != 'infeasible', cutStart:=c(0, head(cutEnd, -1))]
  x$duration[duration == 'infeasible', duration := as.numeric(NA)]
  x$duration$duration = as.numeric(x$duration$duration)
  setkey(x$duration, cut)
  setkey(x$edges, cut)
  x$edges = x$edges[x$duration[!is.na(duration), list(cut, cutStart)]]
  x$edges[, c('start', 'cutStart') := list(start+cutStart, NULL)]

  ## match NAs for message edges
  x$edges[power < 1, power := as.numeric(NA)]
  x
}

.writeILP_prefix = function(prefix){
  origPrefix = sub('_fixedLP', '', prefix)
  eReduced = new.env()
  load(paste('../mergedData', origPrefix, 'Rsave', sep='.'), envir=eReduced)
  
  ranks = unique(eReduced$reduced$assignments$rank)
  ## eRuntimes = new.env()
  ## runtimes =
  ##   load(paste('../', head(eReduced$reduced$assignments, 1)$date,
  ##              '/merged.Rsave',sep=''), envir=eRuntimes)

  ##!@todo eReduced and eRuntimes have enough info to recreate the schedule?
  edges_inv = eReduced$reduced$edges_inv
  setkey(edges_inv, e_uid)
  vertices = eReduced$reduced$vertices
  rSched = eReduced$reduced$schedule
  globals = eReduced$reduced$globals
  rm(eReduced)

###!@todo this should be done in loadAndMergeILP, but requires
###!knowledge of vertex order that lives in merged_.Rsave.
  setkey(vertices, vertex)
  orderedCuts =
    vertices[J(
      resultsFixedLPMerged[[prefix]][[1]]$duration$cut
      )][order(start), vertex]
  resultsILPMerged =
    lapply(resultsILPMerged, lapply, accumulateCutStarts, orderedCuts)
  resultsFixedLPMerged =
    lapply(resultsFixedLPMerged, lapply, accumulateCutStarts, orderedCuts)
  assign('resultsILPMerged', resultsILPMerged, envir=.GlobalEnv)
  assign('resultsFixedLPMerged', resultsFixedLPMerged, envir=.GlobalEnv)
  cols =
    c('src', 's_uid', 'd_uid', 'dest', 'type', 'start', 'weight',
      ## name
      'size',
      ## dest
      ## src
      'tag',
      'power', 'OMP_NUM_THREADS', 'cpuFreq')
  vertexCols = c('vertex', 'label', 'hash', 'reqs')

  .writeILP_prefix_powerLimit = function(pl){
    if(!nrow(pl$edges)){
      cat("powerLimit", pl$duration$powerLimit, ": no edges\n")
      return(NULL)
    }
    setkey(pl$edges, e_uid)

    ## merge sched with reduced edges_inv
    sched = edges_inv[pl$edges]
    save(sched,
         file=
         paste('sched', prefix,
               paste('p', pl$duration$powerLimit[1], 'w', sep=''),
               'Rsave', sep='.'))

    .writePowerTime = function(s, label){
      write.table(powerTime(s),
                  file=
                  paste('powerTime',
                        prefix, label, paste('p', pl$duration$powerLimit[1], 'w', sep=''),
                        'dat', sep='.'),
                  quote=F, sep='\t', row.names=F)
    }
    .writePowerTime(sched, 'all')
    lapply(ranks, function(r)
           .writePowerTime(sched[rank == r], label=sprintf('%06d', r)))

    ##!@todo plot per-rank power vs time, compare power allocation nonuniformity
    
    cols = intersect(cols, names(sched))
    setkey(sched, src)
    schedDest = data.table::copy(sched)
    setkey(schedDest, dest)
    if(nrow(schedDest[J(2)]) != length(ranks)){
      stop("MPI_Finalize anomaly found in LP schedule  for prefix ", prefix, "!\n")
    }


    .writeILP_prefix_powerLimit_rank = function(r){
      ## rank == dest rank for messages
      compEdges = 
        vertices[, vertexCols, with=F][cbind(
                                 sched[type=='comp' & rank==r,
                                       cols, with=F],
                                 d_rank=as.integer(NA),
                                 s_rank=as.integer(NA),
                                 mseq=as.numeric(NA))]
      if(nrow(sched[type == 'message'])){
        messageSendEdges = sched[type=='message' & s_rank==r]
        messageRecvEdges = sched[type=='message' & rank==r]

        messageSendEdges[, d_rank := rank]
        
        setkey(messageRecvEdges, 'o_dest')
        messageRecvEdges[, d_rank := rank]
        
        messageSendEdges = vertices[, vertexCols,
          with=F][messageSendEdges[, c(cols, 's_rank', 'd_rank'), with=F]]
        messageRecvEdges = vertices[, vertexCols,
          with=F][messageRecvEdges[, c(cols, 's_rank', 'd_rank', 'o_dest', 'o_d_uid'), with=F]]
        messageRecvEdges[, src := NULL]
### reduceConfs produces one message edge for each send/recv
### pair, but we want a separate row for both send and recv
        messageEdges =
          .rbindlist(list(messageRecvEdges,
                          cbind(messageSendEdges, o_d_uid=as.numeric(NA))))
        messageEdges[,mseq:=max(s_uid, o_d_uid, na.rm=T),by=vertex]
        messageEdges[, o_d_uid := NULL]
        
        edges =
          .rbindlist(list(compEdges,
                          messageEdges))
        rm(messageEdges)
      } else
        edges = compEdges
      edges[, mseq:=max(mseq, s_uid, na.rm=T), by=list(vertex,type)]
      
      edges = 
        edges[,
              if(.N > 1){
### we should never have more than one comp edge and one message edge
### leaving a vertex
                a = .SD[type=='comp']
                a[, label := as.character(NA)]
                a =
                  rbindlist(list(cbind(.SD[type=='message'], seq=1),
                                 cbind(a, seq=2)))
### hack to handle start times from sender in recv edges.
                a[, start := min(start)]
              } else {
                a=copy(.SD)
                a[,c('label', 'hash', 'reqs'):=as.character(NA)]
                cbind(rbindlist(list(.SD,a)), seq=1)
              },
              by=vertex]

      edges[,mseq:=min(mseq),by=list(vertex)]
      edges = edges[order(mseq, seq)]

      ##!@todo UMT is missing MPI_Finalize; WTF?
      ## handle finalize
      edges =
        .rbindlist(
          list(
            edges,
            vertices[, vertexCols,
                     with=F][cbind(
                       schedDest[dest==2 & rank==r, cols, with=F],
                       d_rank=as.integer(NA), s_rank=as.integer(NA),
                       mseq=max(edges$mseq) + 1, seq=1)]))
      
      edges[, c('seq', 'vertex', 's_uid') := NULL]
      edges[, c('src', 'dest'):=as.integer(NA)]
      edges[type == 'message',
            c('src', 'dest') := list(as.integer(s_rank), as.integer(d_rank))]
      edges[, c('d_rank', 'type') := NULL]
      if(!'size' %in% names(edges)){
        edges[, c('size', 'tag', 'comm') := as.numeric(NA)]
      }
      edges =
        edges[, list(start,
                     duration=weight,
                     name=sapply(strsplit(label, ' '), '[[', 1),
                     size,
                     dest,
                     src,
                     tag,
                     comm='0x0', ##!@todo fix
                     hash,
                     flags=0, ##!@todo fix?
                                        #pkg_w=power,
                                        #pp0_w=0,
                                        #dram_w=0,
                     reqs,##=as.character(NA), ##!@todo fix
                     OMP_NUM_THREADS,
                     cpuFreq
                     )]
      for(col in setdiff(names(edges), c('reqs', 'name', 'comm', 'hash', 'pkg_w'))){
        eCol = edges[[col]]
        eCol[is.na(eCol)] = globals$MPI_UNDEFINED
        edges[[col]] = eCol
      }
      edges[is.na(comm), comm:=0]
      edges[is.na(hash), hash:='0']
                                        #edges[is.na(pkg_w), pkg_w:=0]
      edges[, cpuFreq:=as.integer(cpuFreq)]
      edges[!is.na(name), duration := 0.0]
      edges[, reqs:=sapply(reqs, paste, collapse=',')]
      write.table(edges,
                  ## C code uses %s.%06d.dat
                  file=
                  paste('replay', prefix,
                        paste('p', pl$duration$powerLimit[1], 'w', sep=''),
                        sprintf('%06d', r),
                        'dat', sep='.'),
                  quote=F, sep='\t', row.names=F)
    }
    
    ##debug(.writeILP_prefix_powerLimit_rank)
    lapply(ranks, .writeILP_prefix_powerLimit_rank)
  }

  mclapply(resultsFixedLPMerged[[prefix]], .writeILP_prefix_powerLimit)
  NULL
}

## For each command, for each power limit, write a configuration
## schedule. This involves matching scheduled edges with edges from
## the original schedule, verifying that all edges were scheduled in
## the solution, matching edges with corresponding start vertices,
## writing vertices and edges in start order per rank, etc. We also
## require request IDs and communicator IDs. Perhaps it would be
## easier to load an existing replay schedule and add config options.
writeILPSchedules = function(){
  nnapply(names(resultsFixedLPMerged), .writeILP_prefix)
}

summarizeSchedules = function(){
  napply(
    resultsFixedLPMerged,
    function(results, name){
      prefix = name
      napply(
        results,
        function(plResults, name){
          powerLimit = name
          powerTime = 
            fread(paste('powerTime', prefix, 'all',
                        paste('p',
                              ###!@todo this should agree with .writePowerTime()
                              as.integer(powerLimit),
                              'w', sep=''),
                        'dat', sep='.'))
          ## plot(stepfun(powerTime$start, c(powerTime$power,0)))
          meanPower =
            sum(powerTime[, diff(start) * tail(power, -1)])/tail(powerTime[, start],1)
          plResults$edges[, list(duration=max(start+weight),
                                 meanPower = meanPower,
                                 maxPower = max(powerTime$power))]
        })})
}

if(!interactive()){
##  loadAndMergeLP()
  loadAndMergeILP()
  writeILPSchedules()
}

