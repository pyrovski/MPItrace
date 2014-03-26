#!/usr/bin/env Rscript

source('~/local/bin/pbutils.R')

require('data.table')
require('parallel')
require('igraph')
require('hash')

adjWidth()

debug=T

if(debug){
  options(mc.cores=1)
}
options(datatable.nomatch=0)

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

MPI_req_sources =
  c('MPI_Isend','MPI_Irecv','MPI_Irsend','MPI_Bsend_init',
    'MPI_Rsend_init','MPI_Send_init','MPI_Ssend_init',
    'MPI_Recv_init', 'MPI_Ibsend', 'MPI_Issend')

MPI_Sends =
  c('MPI_Send', 'MPI_Isend', 'MPI_Issend', 'MPI_Bsend', 'MPI_Rsend',
    'MPI_Irsend', 'MPI_Ibsend', 'MPI_Ssend', 'MPI_Rsend_init',
    'MPI_Ssend_init', 'MPI_Bsend_init', 'MPI_Send_init')

MPI_Recvs =
  c('MPI_Recv', 'MPI_Irecv', 'MPI_Recv_init')

## for Merlot, I have the following formula for message latency:
## 6.456e-06 + 2.478e-10 * size
latency = list(merlot = function(size) 6.456e-6 + 2.478e-10 * size)

readGlobal = function(path = '.', filename = "glog.dat"){
  source(filename)
}

readRuntime = function(filename){
  firstLine = strsplit(readLines(filename, n=1),'\t')[[1]]
  colClasses = c(
    start='numeric',
    duration='numeric',
    name='character',
    size='integer',
    dest='integer',
    src='integer',
    tag='integer',
    comm='character',
    hash='character',
    pkg_w='numeric',
    pp0_w='numeric',
    dram_w='numeric',
    reqs='character')
  a = data.table(read.table(filename,h=T,stringsAsFactors=F,colClasses=colClasses))
  a[size == MPI_UNDEFINED,]$size = NA
  a[dest == MPI_UNDEFINED,]$dest = NA
  a[src == MPI_UNDEFINED,]$src = NA
  a[tag == MPI_UNDEFINED]$tag = NA
  a$reqs = strsplit(a$reqs,',')
  a[comm == '(nil)',comm:='']
  return(a)
}

readAll = function(path='.'){
  readGlobal(path)

  ## runtime
  files = sort(list.files(path,'runtime.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  runtimes = lapply(files, readRuntime)
  names(runtimes) = ranks
  runtimes = napply(runtimes, function(x, name){x$rank = as.numeric(name);x})

  ## assignment
  files = sort(list.files(path,'assign.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  assignments = rbindlist(lapply(files, read.table, h=T))
  
  return(list(runtimes=runtimes, assignments=assignments))
}

## single rank only!
.deps = function(x, maxRank){
  ranks = unique(x$rank)
  if(length(ranks) > 1)
    stop('Too many ranks in deps()!')

  rank = ranks
  rm(ranks)

  ##! process ANY_SOURCE and ANY_TAG entries
  ##!@todo preserve ANY_SOURCE and ANY_TAG in new columns?
  ##!@todo resolve ANY_SOURCE and ANY_TAG for Irecv (not a problem in MCB)
  sel = x$src == MPI_ANY_SOURCE | x$tag == MPI_ANY_TAG
  ##!@todo for replay, we need to know both the ANY_SOURCE for
  ##! matching and the actual source to post the PMPI receive.
  if(length(sel))
    x = x[!sel]

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
      newRank = x[s+1,tag]
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
  sources =
    x[name %in% MPI_req_sources, which=T]
  ## only entries with requests
  sources =
    intersect(sources,
              which(sapply(x$reqs, function(x) !all(is.na(x)))))
  
  sinks =
    x[name %in% MPI_req_sinks, which=T]
  ## only entries with requests
  sinks =
    intersect(sinks,
              which(sapply(x$reqs, function(x) !all(is.na(x)))))

  reqs = sort(c(sinks, sources))

  if(length(sinks)){
    sinkReqs = unique(unlist(x[sinks]$reqs))
    
    ## for each request sink, find its source
    ## double-link sinks and sources for later matching of messages
    sinkIndex = 1
    for(req in sinkReqs){
      if(debug && !(sinkIndex %% 100))
        cat('Rank', rank, 'sink request', sinkIndex, 'of', length(sinkReqs), '\n')
      reqFound =
        reqs[sapply(x[reqs]$reqs, function(x) req %in% x)]
###      cat('Rank', rank, 'reqFound:', reqFound, '\n')
      sourced = F
      for(i in reqFound){
        if(x[i]$name %in% MPI_req_sources){
          sourced = T
          last = i
        } else if(sourced){ ## sink
          x$brefs[[i]] = union(unlist(x[i]$brefs), x$uid[last])
          x[last]$fref = x$uid[i]
          sourced = F
        }
      }
      sinkIndex = sinkIndex + 1
    }
  }
  x[,reqs:=NULL]
  cat('Rank', rank, 'done\n')
  return(list(runtimes = x, comms = commTable))
}

##!@todo make sure non-WORLD comm messages are matched
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
    return(list(runtimes=x, comms=commTable))
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
  srDeps = mclapply(1:nrow(mids), function(mid){
    matching = x[mids[mid]][order(uid)]
    canceled =
      x[uid %in% matching[!is.na(fref), fref] & name == 'MPI_Cancel', uid]
### ignore canceled requests and requests not waited for
    matching = matching[!fref %in% canceled & !is.na(fref)]
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
    if(nrow(sends) != nrow(recvs))
      stop('mismatched send-receive:', mids[mid])

    result =
      lapply(1:nrow(sends), function(row) f_noSideEffects(sends[row], recvs[row]))
    result = as.data.frame(do.call(rbind, result))
    names(result) = c('src', 'dest')
    return(result)
  })
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

deps = function(x){
  ranks = sapply(x, function(x) unique(x$rank))
  
  ## if(debug)
  ##   x = lapply(x, .deps, max(ranks))
  ## else
    x = mclapply(x, .deps, maxRank=max(ranks))
  
  ## merge tables
  commTable = rbindlist(lapply(x, '[[', 'comms'))
 
  ##!@todo order by topological sort after all dependecies done
  x = rbindlist(lapply(x, '[[', 'runtimes'))[order(start)]

  ## remap uids
  newUIDs = as.list(1:nrow(x))
  names(newUIDs) = x$uid
  x$uid = unlist(newUIDs[as.character(x$uid)])

  sel = x[!is.na(fref), which=T]
  x$fref[sel] = unlist(newUIDs[as.character(x$fref[sel])])

  sel = x[!is.na(deps), which=T]
  x$deps[sel] = newUIDs[as.character(x$deps[sel])]

  sel = x[!is.na(succ), which=T]
  x$succ[sel] = newUIDs[as.character(x$succ[sel])]

  sel = which(!sapply(x$brefs, is.null))
  ##!@todo speed this up
  x$brefs[sel] =
    lapply(x$brefs[sel], function(x)
           unlist(unname(newUIDs[sapply(x,as.character)])))

  if(nrow(commTable) > 0){
    commTable$source = newUIDs[as.character(commTable$source)]

    sel = which(!is.na(commTable$sink))
    commTable$sink[sel] = newUIDs[as.character(commTable$sink[sel])]
    
    commTable$source = unlist(commTable$source)
    commTable$sink = unlist(commTable$sink)

    ## For now, assume only a single level of derived communicators.
    if(any(commTable$parentComm != MPI_COMM_WORLD))
      cat('Multiple-derived communicators not supported yet.\n')

    ## unify communicators
    commTable$done = F
    cid = 1
    
    commList = list(MPI_COMM_WORLD)
    while(length(commList)){
      comm = commList[[1]]
      commList = tail(commList, -1)
      setkey(commTable, rank, source, done)
      sel = commTable[parentComm == comm & !done, which=T]

      d = commTable[sel,list(source=min(source)),by=rank]
      d$done = F
      setkey(d)
      sel = commTable[d, which=T]
      sel = intersect(sel, commTable[childComm != MPI_COMM_NULL, which=T])

      ## Handle multiple subcomms from the same collective. The
      ## selection represents a single collective call, but the
      ## resulting communicators may be disjoint
      uranks = unique(commTable[sel, ranks])
      commTable[sel, unifiedComm :=
                as.character(cid-1+match(commTable[sel,ranks], uranks))]
      commTable[d, done := T]
      cid = max(as.integer(commTable$unifiedComm), na.rm=T) + 1
      
      ## replace newly unified parentComms in commTable
      setkey(commTable, rank, parentComm)
      commMap =
        unique(commTable[!is.na(unifiedComm), list(rank, childComm, unifiedComm)])
      setnames(commMap, names(commMap), c('rank', 'parentComm', 'unifiedComm'))
      setkey(commMap, rank, parentComm)
      commTable[commMap, parentComm := unifiedComm]

      ## replace childComms in x
      
      ##! this handles multiple uses of the same communicator value in
      ##! different communicator lifetimes, and makes sure replacement
      ##! locations are between source and sink (inclusive)
      d$done = T
      setkey(d)
      setkey(commTable, rank, source, done)
      rowApply(commTable[d], function(row)
               x[rank == row$rank &
                 comm == row$childComm &
                 uid >= row$source &
                 uid <= row$sink,
                 comm := row$unifiedComm]
               )

      ## if more entries in current comm, continue
      d$source = NULL
      d$done = F
      setkey(d)
      setkey(commTable, rank, done)
      if(nrow(commTable[d]))
        commList = c(comm, commList)

      ##!@todo get next comm (must be a child comm by definition
    }
    commTable$done = NULL
  } ## if commTable not empty
  
  ## collectives
  ## init and finalize
  init = which(x$name == 'MPI_Init' | x$name == 'MPI_Init_thread')
  finalize = which(x$name == 'MPI_Finalize')
  vid = 3

  ## g = data.table(name='MPI_Init', succ=as.integer(NA), vertex=1)
  x[init, 'vertex'] = 1
  
  ## g = rbindlist(list(g, list(name='MPI_Finalize', succ=NA, vertex=2)))
  x[finalize, 'vertex'] = 2

  if(debug)
    cat('Collectives\n')
  collectives = intersect(x$name, setdiff(MPI_collectives, c('MPI_Init','MPI_Finalize')))
  
  setkey(x, name, comm, rank)
  for(a_coll in collectives){
    instances = which(x$name == a_coll)
    u = unique(x[instances, list(name, comm)])
    for(rank in ranks){
      urank = c(u, rank=rank)
      x[urank, vertex:=vid+(0:(nrow(.SD)-1))]
      vidInc = nrow(x[urank])
    }
    vid = vid + vidInc
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

### vertices and edges with attributes.
tableToGraph = function(x, assignments){
  ## replace computation vertices with weighted edges
  setkey(x, uid)
  uids = x[is.na(name)]$uid
  edges = do.call(rbind, lapply(uids,
    function(u) {
      e=x[J(u)];
      data.frame(src=x[J(unlist(e$deps))]$vertex,
                 dest=x[J(unlist(e$succ))]$vertex,
                 weight=e$duration)
    }))
  ##predMap = hash(uids, edges$src)
  ##succMap = hash(uids, edges$dest)

  x = x[!is.na(name)]
  
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
    f2 = function(d) data.frame(src=x[J(d)]$vertex, dest=x[J(u)]$vertex)
    rbindlist(lapply(x[J(u)]$deps, f2))
  }

  ## add communication edges
  edges = rbind(edges, cbind(rbindlist(lapply(sel, f)), weight=0))
    
  g = graph.data.frame(edges, vertices=vertices)
  
  return(g)
}

run = function(path='.'){
  a = readAll(path)
  b = messageDeps(deps(a$runtimes))
  return(b)
}
