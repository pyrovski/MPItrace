#!/usr/bin/env Rscript

##!@todo use hash tables for multi-element columns? Apparently,
## data.table in newer versions of R refuses to set keys on a data
## table that contains list columns. There's a test with
## isVectorAtomic() in reorder.c that fails on lists, even though the
## documentations says lists are supported.

source('~/local/bin/pbutils.R')

require('data.table')
require('parallel')
require('igraph')
require('hash')

adjWidth()

debug=T

if(debug){
##  options(mc.cores=1)
}
options(datatable.nomatch=0)

colClasses = c(
  start='numeric',
  duration='numeric',
  name='character',
  size='integer',
  dest='integer',
  src='character', ## converted later to integer
  tag='character', ## converted later to integer
  comm='character',
  hash='character',
  pkg_w='numeric',
  pp0_w='numeric',
  dram_w='numeric',
  reqs='character' ## converted later to character list
  )

MPI_Comm_sources =
  c('MPI_Comm_dup', 'MPI_Comm_create', 'MPI_Comm_split')

MPI_Comm_sinks =
  c('MPI_Comm_free')

MPI_collectives = c(
  'MPI_Init',
  'MPI_Finalize',
  'MPI_Barrier',
  'MPI_Bcast',
  'MPI_Scatter',
  'MPI_Scatterv',
  'MPI_Gather',
  'MPI_Gatherv',
  'MPI_Allgather',
  'MPI_Allgatherv',
  'MPI_Reduce',
  'MPI_Allreduce',
  'MPI_Reduce_scatter',
  'MPI_Alltoall',
  'MPI_Alltoallv',
  'MPI_Scan',
  'MPI_File_write_at_all',
  MPI_Comm_sources)

MPI_req_sinks =
  c('MPI_Wait', 'MPI_Waitall', 'MPI_Waitsome',
    'MPI_Test', 'MPI_Testall', 'MPI_Testsome',
    'MPI_Request_free', 'MPI_Cancel')

MPI_Isends =
  c('MPI_Isend','MPI_Irsend','MPI_Bsend_init',
    'MPI_Rsend_init','MPI_Send_init','MPI_Ssend_init',
    'MPI_Ibsend','MPI_Issend')

MPI_Irecvs =
  c('MPI_Irecv','MPI_Recv_init')

MPI_req_sources =
  c(MPI_Isends, MPI_Irecvs)

MPI_Sends =
  c('MPI_Send', 'MPI_Bsend', 'MPI_Rsend',
    'MPI_Ssend', MPI_Isends)

MPI_Recvs =
  c('MPI_Recv', MPI_Irecvs)

## for Merlot, I have the following formula for message latency:
## 6.456e-06 + 2.478e-10 * size
latency = list(merlot = function(size) 6.456e-6 + 2.478e-10 * size)

readGlobal = function(path = '.', filename = "glog.dat"){
  source(file.path(path, filename))
}

readRuntime = function(filename, path='.'){
  filename = file.path(path, filename)
  firstLine = strsplit(readLines(filename, n=1),'\t')[[1]]
  a = data.table(read.table(filename,h=T,stringsAsFactors=F,colClasses=colClasses))
  a[size == MPI_UNDEFINED,]$size = NA
  a[dest == MPI_UNDEFINED,]$dest = NA

  a$src = lapply(strsplit(a$src, ','), as.integer)
  a$src = lapply(a$src, function(x) {x[x == MPI_UNDEFINED] = NA; x})
  
  a$tag = lapply(strsplit(a$tag, ','), as.integer)
  a$tag = lapply(a$tag, function(x) {x[x == MPI_UNDEFINED] = NA; x})

  a$reqs = strsplit(a$reqs,',')
  ##a[comm == '(nil)',comm:='']
  return(a)
}

readAll = function(path='.'){
  ##!@todo for comparing multiple runs, the globals may be different.
  readGlobal(path)

  ## runtime
  files = sort(list.files(path,'runtime.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  runtimes = lapply(files, readRuntime, path=path)
  names(runtimes) = ranks
  runtimes = napply(runtimes, function(x, name){x$rank = as.numeric(name);x})

  ## assignment
  files = sort(list.files(path,'assign.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  assignments =
    rbindlist(lapply(files, function(file)
                     read.table(file.path(path, file), h=T)))
  
  return(list(runtimes=runtimes, assignments=assignments))
}

## single rank only!

## Writing replays is important because the resulting traces will have
## resolved MPI_ANY_TAG and MPI_ANY_SOURCE.
.deps = function(x, maxRank, writeReplay=T, path='.'){
  ranks = unique(x$rank)
  if(length(ranks) > 1)
    stop('Too many ranks in deps()!')

  rank = ranks
  rm(ranks)

  ##! process ANY_SOURCE and ANY_TAG entries
  ##!@todo preserve ANY_SOURCE and ANY_TAG in new columns?
  ##!@todo resolve ANY_SOURCE and ANY_TAG for Irecv (not a problem in MCB)
  ##sel = (x$src == MPI_ANY_SOURCE | x$tag == MPI_ANY_TAG) & x$name == 'MPI_Recv'
  ##!@todo for replay, we need to know both the ANY_SOURCE for
  ##! matching and the actual source to post the PMPI receive?
  ## if(length(which(sel))){
  ##   if(debug){
  ##     cat('Rank', rank, 'deleting', length(which(sel)), 'duplicate ANY_SOURCE and ANY_TAG entries:\n')
  ##     print(x[sel])
  ##   }
  ##   x = x[!sel]
  ## }

  x$vertex = as.numeric(NA)

  # numeric for now, list later
  ## deps: predecessors
  x$deps = as.numeric(NA)

  ## succ: successors (only maintained per rank)
  x$succ = as.numeric(NA)
  
  ## brefs: request sources
  x$brefs = vector('list', nrow(x))

  ## fref: request sink
  x$fref = as.numeric(NA)

  ## uid: unique identifier
  x$uid = (1:nrow(x))+(rank+1)/(maxRank+1)

  ## clean up duplicate communicator entries for comm creation
  ## tag: new rank
  ## msgSize: new size
  ## comm: new comm
  ## if not participating, new comm = MPI_COMM_NULL (will not free)
  if(any(x$name %in% MPI_Comm_sources)){
    ## remove duplicate communicator creation rows, add to a separate table

    ## old comm, new comm, new size, new rank, uid of sink

    sel =
      x[name %in% MPI_Comm_sources & is.na(size) & comm != MPI_COMM_NULL,
        which=T]

    f = function(s) {
      parentComm = x[s,comm]
      childComm = x[s+1,comm]
      newSize = x[s+1,size]
      newRank = unlist(x[s+1,tag])
      source = x[s, uid]
      if(childComm == MPI_COMM_NULL)
        sink = as.numeric(NA)
      else {
        sink = head(x[uid > x[s, uid] &
          name %in% MPI_Comm_sinks &
          comm == x[s+1,comm],
          uid],
          1)
        if(length(sink) < 1)
          sink = x[name == 'MPI_Finalize' & rank == x[s, rank], uid]
      }
      return(data.frame(parentComm = parentComm,
                        childComm = childComm,
                        newSize = newSize,
                        newRank = newRank,
                        source = source,
                        sink = sink,
                        stringsAsFactors=F))
    }
    commTable = rbindlist(lapply(sel, f))
    commTable$rank = rank
    commTable$ranks =
      lapply(commTable$source, function(r)
             as.integer(unlist(x[1 + x[uid == r, which=T], reqs])))

    sel = x[!name %in% MPI_Comm_sources | is.na(size) & comm != MPI_COMM_NULL, which=T]
    x = x[sel]
  } else
    commTable = NULL

  ## sequential dependencies
  rows = 1:nrow(x)
  x$deps[tail(rows, -1)] = x$uid[head(rows, -1)]
  x$succ[head(rows, -1)] = x$uid[tail(rows, -1)]

  ## get backrefs for blocking receives and waits,
  ## nonblocking successful tests, and request frees
  if(debug)
    cat('Rank', rank, 'requests\n')

  ## only entries with requests
  setkey(x, uid)
  reqUIDs = x[!is.na(reqs), uid]
  sources = intersect(x[name %in% MPI_req_sources, uid], reqUIDs)
  sinks = intersect(x[name %in% MPI_req_sinks, uid], reqUIDs)
  
  if(length(sinks)){
    sinkReqs = unique(unlist(x[J(sinks)]$reqs))

    ##!@todo speed this up. I need to find source->sink pairs
    ##!quickly. We could use mclapply here if we merged the brefs
    ##!lists after the fact, and set fref, tag, and src.
    
    ## for each request sink, find its source
    ## double-link sinks and sources for later matching of messages
    sinkCount = 1
    for(req in sinkReqs){
      if(debug && !(sinkCount %% 100))
        cat('Rank', rank, 'sink request', sinkCount, 'of', length(sinkReqs), '\n')
      ## uids
      reqFound =
        reqUIDs[sapply(x[J(reqUIDs)]$reqs, function(x) req %in% x)]
      sourced = F
      for(i in reqFound){
        sinkIndex = x[J(i), which=T]
        if(x[sinkIndex]$name %in% MPI_req_sources){
          sourced = T
          last = i ## uid
        } else if(sourced){ ## sink
          sourceIndex = x[J(last), which=T]
          x[sinkIndex, list(brefs = as.list(union(unlist(x[sinkIndex]$brefs), last)))]
          x[sourceIndex, fref := i]
          if(x[sourceIndex]$tag == MPI_ANY_TAG ||
             x[sourceIndex]$src == MPI_ANY_SOURCE){
            matchIndex = which(x[sinkIndex]$reqs == req)
            resolvedTag = x[sinkIndex]$tag[matchIndex]
            resolvedSrc = x[sinkIndex]$src[matchIndex]
            x[sourceIndex, list(tag=resolvedTag, src=resolvedSrc)]
            cat('resolved tag and source for req', req, '\n')
          }
          sourced = F
        }
      }
      sinkCount = sinkCount + 1
    }
  }
  if(!debug)
    x[,reqs:=NULL]
  ## remove src and tag from req sinks, convert back from list
  x$src = unlist(rowApply(x, function(row){
    if(row$name %in% MPI_req_sinks)
      as.integer(NA)
    else
      row$src
  }))
  x$tag = unlist(rowApply(x, function(row){
    if(row$name %in% MPI_req_sinks)
      as.integer(NA)
    else
      row$tag
    }))
  if(writeReplay){
    cat('Rank', rank, 'writing replay\n')
    ##!@todo the replay fscanf is going to trip on NAs
    y = data.table::copy(x[, names(colClasses), with=F][order(start)])
    y$reqs = unlist(lapply(y$reqs, paste, collapse=','))
    write.table(y,
                file=
                file.path(path,
                          paste('replay', sprintf('%06d', rank), 'dat', sep='.')),
                quote=F, sep='\t', row.names=F, na=as.character(MPI_UNDEFINED))
    rm(y)
  }
  cat('Rank', rank, 'done\n')
  return(list(runtimes = x, comms = commTable))
}

preDeps = function(x, ...){
  ranks = sapply(x, function(x) unique(x$rank))
  
  x = mclapply(x, .deps, maxRank=max(ranks), ...)
  return(x)
}

deps = function(x){
  if(debug)
    cat('Merging tables and remapping UIDs\n')
  ## merge tables
  commTable = rbindlist(lapply(x, '[[', 'comms'))
 
  x = rbindlist(lapply(x, '[[', 'runtimes'))[order(start)]
  ranks = sort(unique(x$rank))

  ## remap uids
  newUIDs = as.list(1:nrow(x))
  names(newUIDs) = x$uid
  x$uid = unlist(newUIDs[as.character(x$uid)])

  sel = x[!is.na(fref), which=T]
  x$fref[sel] = unlist(newUIDs[as.character(x$fref[sel])])

  sel = x[!is.na(deps), which=T]
  x$deps[sel] = newUIDs[as.character(x$deps[sel])]

  sel = x[!is.na(succ), which=T]
  x$succ[sel] = unlist(newUIDs[as.character(x$succ[sel])])

  sel = which(!sapply(x$brefs, is.null))
  ##!@todo speed this up
  x$brefs[sel] =
    lapply(x$brefs[sel], function(x)
           unlist(unname(newUIDs[sapply(x,as.character)])))

  if(debug)
    cat('Unifying communicators\n')
  if(nrow(commTable) > 0){
    commTable$source = newUIDs[as.character(commTable$source)]

    sel = which(!is.na(commTable$sink))
    commTable$sink[sel] = newUIDs[as.character(commTable$sink[sel])]
    
    commTable$source = unlist(commTable$source)
    commTable$sink = unlist(commTable$sink)

    ## unify communicators
    commTable$done = F
    cid = 1
    
    commList = list(MPI_COMM_WORLD)
    while(length(commList)){
      comm = commList[[1]]
      commList = tail(commList, -1)
      setkey(commTable, rank, source, done)
      sel = commTable[parentComm == comm & !done, which=T]

      ##!@todo this doesn't handle terminal child comms yet
      d = commTable[sel,list(source=min(source)),by=rank]
      d$done = F
      setkey(d)
      sel = commTable[d, which=T]
      sel = intersect(sel, commTable[childComm != MPI_COMM_NULL, which=T])

      ## Handle multiple subcomms from the same collective. The
      ## selection represents a single collective call, but the
      ## resulting communicators may be disjoint
      uranks = unique(commTable[sel, ranks])
      oldUnified = Filter(Negate(is.na), unique(commTable$unifiedComm))
      commTable[sel, unifiedComm :=
                as.character(cid-1+match(commTable[sel,ranks], uranks))]
      cid = max(as.integer(commTable$unifiedComm), na.rm=T) + 1
      newUnified =
        setdiff(Filter(Negate(is.na), unique(commTable$unifiedComm)), oldUnified)
      
      ## replace newly unified parentComms in commTable
      commMap =
        unique(commTable[!is.na(unifiedComm) & !done, list(rank, childComm, unifiedComm)])
      setnames(commMap, names(commMap), c('rank', 'comm', 'unifiedComm'))
      setkey(commMap, rank, comm)

      ##!@todo use d for selection
      ## childSel = commTable[d, which=T]
      ## if(length(childSel)){
      ##   setkey(commTable, rank, childComm)
      ##   commTable[childSel]$childComm =
      ##     commMap[commTable[childSel, list(rank, comm=childComm)]]$unifiedComm
      ## }

      f = function(row){
        sel =
          commTable[parentComm == row$childComm &
                    source >= row$source &
                    source <= row$sink, which=T]
        if(length(sel))
          commTable[sel, parentComm := row$unifiedComm]
      }
      setkeyv(commTable, key(d))
      childSel = commTable[d, which=T]

      rowApply(commTable[d], f)
      
      ## parentSel = c();
      ## if(length(parentSel)){
      ##   setkey(commTable, rank, parentComm)
      ##   commTable[parentSel]$parentComm =
      ##     commMap[commTable[parentSel, list(rank, comm=parentComm)]]$unifiedComm
      ## }

      setkeyv(commTable, key(d))
      commTable[d, done := T]

      ## replace childComms in x
      
      ##! this handles multiple uses of the same communicator value in
      ##! different communicator lifetimes, and makes sure replacement
      ##! locations are between source and sink (inclusive)
      d$done = T
      setkey(d)
      setkey(commTable, rank, source, done)
      translate = commTable[d]
      setkey(commMap, rank, unifiedComm)
      setkey(translate, rank, childComm)
      ##!@todo test
      translate[commMap, childComm := comm]
      f = function(row)
        x[uid >= row$source & uid <= row$sink & comm  == row$childComm,
          comm := row$unifiedComm]
      rowApply(translate, f)
      rm(translate)
    
      ## if more entries in current comm, continue
      if(nrow(commTable[parentComm == comm & !done]))
        commList = c(comm, commList)
      
      ##!@todo get next comm (must be a child comm by definition)
      derivedParents = intersect(commTable$parentComm, newUnified)
      commList = c(commList, derivedParents)
    }
    commTable$done = NULL
  } ## if commTable not empty
  
  if(debug)
    cat('Unifying collectives\n')
  ## collectives
  ## init and finalize
  vid = 3

  x[name %in% c('MPI_Init', 'MPI_Init_thread'), 'vertex'] = 1
  x[name == 'MPI_Finalize', 'vertex'] = 2

  if(debug)
    cat('Collectives\n')
  collectives = intersect(x$name, setdiff(MPI_collectives, c('MPI_Init','MPI_Finalize')))

  setkey(x, name, comm, rank)
  for(a_coll in collectives){
    instances = which(x$name == a_coll)
    u = unique(x[instances, list(name, comm)])
    setkey(u)
    ## apply to all ranks in a comm at once
    for(i in 1:nrow(u)){
      row = u[i]
      x[row, vertex:=vid+(0:(nrow(.SD)-1)), by=rank]
      vid = max(x[row]$vertex) + 1
    }
  }

  ## add vertices for intra-rank stuff
  x[is.na(vertex), vertex:=vid+(0:(nrow(.SD)-1))]
  vid = max(x$vertex, na.rm=T) + 1

### commTable ###
### parentComm: parent communicator. Mostly useful for collective semantics
### childComm: child communicator
### newSize: size of child communicator
### newRank: resulting rank in child communicator
### source: uid of originating MPI call
### sink: end of communicator validity; MPI_Comm_free() or MPI_Finalize()
### rank: rank of caller in MPI_COMM_WORLD
### ranks: list of MPI_COMM_WORLD ranks participating in new communicator
  return(list(runtimes = x, comms = commTable))
}

messageDeps = function(x){
  ## match inter-rank messages. Messages (excluding MPI_ANY_SOURCE and
  ## MPI_ANY_TAG) are uniquely identified by the following: source,
  ## destination, size, tag, communicator, and order.  For
  ## MPI_ANY_SOURCE and MPI_ANY_TAG, we record the tag and source at
  ## run time.

  commTable = x$comms
  x = x$runtimes
  
  mids = unique(x[name %in% c(MPI_Sends, MPI_Recvs), list(src, dest, size, tag, comm)])
  mids = mids[complete.cases(mids)]

  if(nrow(mids) < 1){
    cat('no messages\n')
    return(x)
  }

  f_noSideEffects = function(s, r){
    if(!is.na(r$fref))
      dest = r$fref
    else
      dest = r$uid
    return(list(src=s$uid, dest=dest))
  }

  setkey(x, src, dest, size, tag, comm)
  ## rbindlist makes shitty data tables?
  f = function(mid){
    matching = x[mids[mid]][order(uid)]
    canceled =
      x[uid %in% matching[!is.na(fref), fref] & name == 'MPI_Cancel', uid]
### ignore canceled requests and requests not waited for
    matching = matching[!fref %in% canceled & !is.na(fref) | !name %in% MPI_req_sources]
    if(nrow(matching) < 1){
      cat('Message ID', mid, 'of', nrow(mids), ':', 'no messages\n')
      return(data.frame())
    }
    
    sends = matching[name %in% MPI_Sends]
    recvs = matching[name %in% MPI_Recvs]
    if(debug)
      cat('Message ID', mid, 'of', nrow(mids), ':', nrow(sends), 'messages\n')
    ## if a receive has a request, the send dependency leads to the
    ## matching wait, test, or free
    if(nrow(sends) != nrow(recvs)){
      errMsg = paste('mismatched send-receive:', paste(mids[mid], collapse='\t'))
      stop(errMsg)
    }

    result =
      lapply(1:nrow(sends), function(row) f_noSideEffects(sends[row], recvs[row]))
    result = as.data.frame(do.call(rbind, result))
    names(result) = c('src', 'dest')
    return(result)
  }
  srDeps = mclapply(1:nrow(mids), f)
  ## this is slower than rbindlist, but doesn't segfault
  srDeps = do.call(rbind, srDeps)
  srDeps = as.data.table(lapply(srDeps, unlist))
  
  if(debug)
    cat('Done finding message matches\n')

  setkey(x, uid)
  ## find indices with multiple references
  setkey(srDeps, dest)
  destRLE = rle(srDeps$dest)
  multDests = destRLE$values[destRLE$lengths > 1]

  ## complete the indices with single references
  sel = setdiff(srDeps$dest, multDests)
  if(length(sel)){
    d = x$deps
    d[x[J(sel), which=T]] = mapply(c, d[x[J(sel), which=T]], srDeps[J(sel)]$src, SIMPLIFY=F)
    x$deps = d
    rm(d)
  }

  ## complete the indices with multiple references
  if(length(multDests)){
    d = x$deps
    sel = x[J(multDests), which=T]
    d[sel] =
      mapply(c, d[sel],
             lapply(multDests,
                    function(d) srDeps[J(d)]$src), SIMPLIFY=F)
    x$deps = d
    rm(d)
  }
    
  if(debug)
    cat('Done matching messages\n')

  return(x)
}

### vertices and edges with attributes.
tableToGraph = function(x, assignments){
  if(debug)
    cat('Computation edges\n')
  ## replace computation vertices with weighted edges
  setkey(x, uid)
  y = data.table::copy(x)
  setkey(y, vertex)
  uids = x[is.na(name), uid]
  edges = rbindlist(mclapply(uids,
    function(u) {
      e=x[J(u)];
      result =
        data.frame(src=x[J(unlist(e$deps))]$vertex,
                   dest=x[J(e$succ)]$vertex,
                   weight=e$duration)
      if(any(is.na(y[J(result$src)]$name))){
        cat('comp->comp!\n')
        print(e)
        print(y[J(c(result$src, result$dest))])
        print(result)
      }
      return(result)
    }))
  ##predMap = hash(uids, edges$src)
  ##succMap = hash(uids, edges$dest)

  x = x[!is.na(name)]

  if(debug)
    cat('Deleting computation predecessors\n')
  ## delete dependencies on computation vertices
  x$deps = mclapply(x$deps, function(e) x[J(e)]$uid)

  ##!@todo decide how to handle collectives (decompose, single vertex, etc.)

  ##!@todo add edge weights to communication edges
  
  ## For the vertex frame, column 1 is the vertex name.
  ##vertices = x[, list(name, size, dest, src, tag, comm, hash, vertex)]
  vertices = x[, list(name, duration, rank, vertex)]
  vertices$name = paste(vertices$name, vertices$rank)
  setnames(vertices, names(vertices), c('label', 'duration', 'rank', 'vertex'))
  vertices$rank = NULL
  ##vertices = x[,list(vertex)]
  setcolorder(vertices, c('vertex', setdiff(names(vertices), 'vertex')))
  setkey(vertices, vertex)
  vertices = unique(vertices)

  ## For the edge list frame, column 1 is the source and column 2 is the dest
  setkey(x, uid)
  sel = x[sapply(x$deps, length) > 0, uid]

  f = function(u){
    xu = x[J(u)]
    uDest = xu$vertex
    f2 = function(d) data.frame(src=x[J(d)]$vertex, dest=uDest)
    rbindlist(lapply(xu$deps, f2))
  }

  if(debug)
    cat('Communication edges\n')
  ## add communication edges
  edges = rbind(edges, cbind(rbindlist(lapply(sel, f)), weight=0))
    
  if(debug)
    cat('Graph object\n')
  g = graph.data.frame(edges, vertices=vertices)
  
  return(g)
}

tableToMarkov = function(x, rank=0){
  .rank = rank
  setkey(x, rank)
  x = x[rank == .rank]
  setkey(x, hash)
  y = data.table::copy(x)
  setkey(y, uid)
  hashes = unique(x[!is.na(name) & name != 'MPI_Finalize']$hash)
  vertices =  x[!is.na(name), list(label=unique(name)), by=hash]
  edges =
    rbindlist(
      mclapply(hashes, function(h){
      ## successor hashes
      successors = y[J(y[J(x[J(h)]$succ)]$succ)]$hash
      if(length(successors)){
        successors = table(successors)
        successors = 10*successors/sum(successors)
        rbindlist(
          lapply(1:length(successors), function(i)
                 data.frame(src=h, dest=names(successors)[i],
                            penwidth=successors[i], stringsAsFactors=F)))
      } else
        data.frame(src=character(0), dest=character(0), weight=numeric(0))
    })
      )
  g = graph.data.frame(edges, vertices=vertices)
  write.graph(g, file='markov.dot', format='dot')
  return(g)
}

shortStats = function(x, thresh=.001){
  shortTaskRatio = nrow(x[duration < thresh])/nrow(x)
  shortTimeRatio = sum(x[duration < thresh, duration])/sum(x$duration)
  cat(shortTaskRatio * 100, '% of tasks <', thresh, 's\n')
  cat(shortTimeRatio * 100, '% of time in short tasks\n')
}

run = function(path='.'){
  a = readAll(path)
  assignments = a$assignments
  b = preDeps(a$runtimes)
  rm(a)
  b = deps(b)
  comms = b$comms
  b2 = messageDeps(b)
  rm(b)
  g = tableToGraph(b2)
  
  return(list(runtimes = b2,
              graph = g,
              assignments = assignments,
              comms = comms))
}
