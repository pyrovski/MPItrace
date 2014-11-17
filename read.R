#!/usr/bin/env Rscript

##!@todo use hash tables for multi-element columns? Apparently,
## data.table in newer versions of R refuses to set keys on a data
## table that contains list columns. There's a test with
## isVectorAtomic() in reorder.c that fails on lists, even though the
## documentations says lists are supported. This is not present in the
## most recent code (1.93), though.

##!@todo what happens to spin rows?

source('~/local/bin/pbutils.R')

require('data.table')
require('parallel')
require('igraph')

adjWidth()

debug=T

#options(mc.cores=16)
options(datatable.nomatch=0)

#!@todo get this from a file
flagBits = list(omp=1, spin=2, newComm=0x1000) ## bit masks, not indices
minDuration = .0000001

new_MPI_COMM_WORLD = '-1'
new_MPI_COMM_NULL = '-2'
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
  flags='numeric',
  pkg_w='numeric',
  pp0_w='numeric',
  dram_w='numeric',
  reqs='character' ## converted later to character list
  )

MPI_Inits = c('MPI_Init', 'MPI_Init_thread')

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

##!@todo handle MPI_Start, MPI_*_init
MPI_Isends =
  c('MPI_Isend','MPI_Irsend',
    'MPI_Ibsend','MPI_Issend')

MPI_Irecvs =
  c('MPI_Irecv')

##! update shim_pre_1 when changing req_sources
MPI_Req_sources =
  c(MPI_Isends, MPI_Irecvs)

MPI_Send_inits = 
  c('MPI_Bsend_init',
    'MPI_Rsend_init','MPI_Send_init','MPI_Ssend_init')
MPI_Recv_inits = c('MPI_Recv_init')
MPI_Req_inits = c(MPI_Send_inits, MPI_Recv_inits)

## these calls don't generate or terminate the request(s), but
## initiate the communication; they act as both sinks (from _init())
## and sources
MPI_Req_starts =
  c('MPI_Start', 'MPI_Startall')

MPI_Req_sinks =
  c('MPI_Wait', 'MPI_Waitall', 'MPI_Waitsome',
    'MPI_Test', 'MPI_Testall', 'MPI_Testsome',
    'MPI_Cancel')

MPI_Req_frees = 'MPI_Request_free'

MPI_Sends =
  c('MPI_Send', 'MPI_Bsend', 'MPI_Rsend',
    'MPI_Ssend', MPI_Isends)

MPI_Recvs =
  c('MPI_Recv', MPI_Irecvs)

latency =
  list(merlot = function(size) 6.456e-6 + 2.478e-10 * size,
       cab = function(size){
###!@todo this difference could be due to a blocking threshold
###!@todo vectorize
         if(size < 4000)
           1.608e-05
         else
           4.453e-05 + 3.109e-10 * size
       })

selfLatency =
  list(merlot = latency[['merlot']],
       cab = function(size){
         if(size < 4000)
           1.414e-06
         else
           5.768e-05 + 1.722e-10 * size
       })

#e = new.env()
#load('acceptablePowerModel_conf_only.Rsave', envir=e)
#E5_2670_power_conf_only = get('m3', envir=e)
#rm(e)

cpuMap = list(cab = 'E5_2670', merlot = 'E5_2670')

## For two sockets. Eventually, I would like to use one model per
## socket.
activePower =
  list(
    E5_2670 =
    function(threads,
             cpu_freq = 2600000, ## kHz
             mem_freq = 1600, ## MHz
             SMT = F,
             a_L3_access_rate,
             C0_ratio
             ){
      if(SMT)
        cores = ceiling(threads/2)
      else
        cores = threads
      
      mcp = cpu_freq * mem_freq * 1000
      mctp = mcp * threads
      mcr = cpu_freq / (mem_freq * 1000)
      mctr = mcr / threads
      mcptr = mcp / threads      
      mcrtp = mcr * threads

      ## conf only
      x = data.frame(mem_freq, cores, cpu_freq, mctr, mcptr, SMT, mcp, mcrtp)
      return(predict(E5_2670_power_conf_only, x))
    })
## per socket
idlePower = list(E5_2670 = 19.1)

readGlobal = function(path = '.', filename = "glog.dat"){
  source(file.path(path, filename))
  e = new.env()
  source(file.path(path, filename), local=e)
  assign('globals', as.list(e), envir=.GlobalEnv)
}

readRuntime = function(filename, maxRank, path='.'){
  filename = file.path(path, filename)

###!@todo if the disk quota is exceeded, the runtime files will exist
###!but be empty.
  firstLine = strsplit(readLines(filename, n=1),'\t')[[1]]
  a = data.table(read.table(filename,h=T,stringsAsFactors=F,colClasses=colClasses))
  a[size == MPI_UNDEFINED,]$size = NA
  a[dest == MPI_UNDEFINED,]$dest = NA

  f = function(x) {x[x == MPI_UNDEFINED] = NA; x}
  a$src = lapply(strsplit(a$src, ','), as.integer)
  a$src = lapply(a$src, f)
  
  a$tag = lapply(strsplit(a$tag, ','), as.integer)
  a$tag = lapply(a$tag, f)

  rank = as.integer(rev(strsplit(filename, '[.]')[[1]])[2])
  a$uid = (1:nrow(a))+(rank)/(maxRank+1)
  
  # handle new comms
  sel = bitwAnd(a$flags, flagBits$newComm) > 0
  if(any(sel))
    newComms = which(sel)
  else
    newComms = as.integer(NA)
  
  splitReqs = strsplit(a$reqs[!sel],',')
  uniqueReqs = Filter(Negate(is.na), unique(unlist(splitReqs)))
  uniqueReqs = setdiff(uniqueReqs, MPI_REQUEST_NULL)
  uidsByReq =
    nnapply(uniqueReqs, function(req){
      sel = grep(paste(req, '(,|$)', sep=''), a$reqs)
      a[sel, uid]
    })
  a$reqs[!sel] = splitReqs
  ##a[comm == '(nil)',comm:='']
  return(list(runtime=a, uidsByReq=uidsByReq, newComms=newComms))
}

readAll = function(path='.'){
  ##!@todo for comparing multiple runs, the globals may be different.
  readGlobal(path)

  ## runtime
  files = sort(list.files(path,'runtime.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  runtimes = mclapply(files, readRuntime, maxRank=max(ranks), path=path)
  names(runtimes) = ranks
  uidsByReq = lapply(runtimes, '[[', 'uidsByReq')
  newComms = lapply(runtimes, '[[', 'newComms')
  runtimes = lapply(runtimes, '[[', 'runtime')
  runtimes =
    napply(runtimes, function(x, name){x$rank = as.numeric(name);x})

  ## assignment
  files = sort(list.files(path,'assign.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  assignments =
    .rbindlist(lapply(files, function(file)
                     read.table(file.path(path, file), h=T,
                                stringsAsFactors=F)))
  return(list(runtimes=runtimes, assignments=assignments, uidsByReq=uidsByReq,
              newComms=newComms)) }

## single rank only!

## Writing replays is important because the resulting traces will have
## resolved MPI_ANY_TAG and MPI_ANY_SOURCE.
.deps = function(x, rank, uidsByReq, newComms, writeReplay=F, path='.'){
  ranks = unique(x$rank)
  if(length(ranks) > 1)
    stop('Too many ranks in .deps()!')

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

  ## succ: successors (only maintained per rank)
  x$succ = as.numeric(NA)
  
  x[duration == 0, duration := minDuration]
  
  x[comm == MPI_COMM_WORLD]$comm = new_MPI_COMM_WORLD
  x[comm == MPI_COMM_NULL]$comm = new_MPI_COMM_NULL
  
  ## clean up duplicate communicator entries for comm creation
  ## tag: new rank
  ## msgSize: new size
  ## comm: new comm
  ## if not participating, new comm = MPI_COMM_NULL (will not free)
  if(any(x$name %in% MPI_Comm_sources)){
    ## remove duplicate communicator creation rows, add to a separate table

    ## old comm, new comm, new size, new rank, uid of sink
    sel =
      x[name %in% MPI_Comm_sources & is.na(size) & comm != new_MPI_COMM_NULL,
        which=T]

    f = function(s) {
      parentComm = x[s,comm]
      childComm = x[s+1,comm]
      newSize = x[s+1,size]
      newRank = unlist(x[s+1,tag])
      source = x[s, uid]
      if(childComm == new_MPI_COMM_NULL)
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
    commTable = .rbindlist(lapply(sel, f))
    commTable$rank = rank
    commTable$ranks = ##!@todo this is broken.
      lapply(x[commTable$source+1, reqs],
             function(r) as.integer(strsplit(r, ',')[[1]]))

    sel = x[!name %in% MPI_Comm_sources | is.na(size) & comm != new_MPI_COMM_NULL, which=T]
    x = x[sel]
    rm(sel)
  } else
    commTable = NULL

  ## sequential dependencies
  rows = 1:nrow(x)
  ##x$deps[tail(rows, -1)] = x$uid[head(rows, -1)]
  x$succ[head(rows, -1)] = x$uid[tail(rows, -1)]

  ## get backrefs for blocking receives and waits,
  ## nonblocking successful tests, and request frees
  if(debug)
    cat('Rank', rank, 'requests\n')

  ## only entries with requests
  setkey(x, uid)
  reqUIDs = x[!is.na(reqs), uid]
  requests = names(uidsByReq)
  
  if(length(requests)){
    ## brefs: request sources
    ## fref: request sinks
    brefs = vector('list', length=nrow(x))
    frefs = vector('list', length=nrow(x))

    ## for each request sink, find its source

    ## double-link sinks and sources for later matching of messages

###!a request can have multiple 'sinks' in its lifetime; waiting for
###!finished requests is allowed, and persistent requests can be
###!restarted.
    reqStates = c('inactive', 'init', 'active', 'null')

    reqTotals = sapply(uidsByReq, length)
    reqCount = sum(reqTotals)
    reqCumsum = cumsum(reqTotals)
    reqCumsum = reqCumsum - reqCumsum[1]
    reqProgress = function(req, offset) (unname(reqCumsum[req]) + offset)/reqCount
###!@todo this could be parallelized if the brefs and frefs lists were
###!merged afterward
    cat('Tracking request states:\n')
    reqWarning = function()
      warning('uid ', x[reqIndex, uid], ': invalid transition: ',
              as.character(reqState), ', ', x[reqIndex, name])
    for(req in requests){
      offset = 0
      progress = reqProgress(req, offset) * 100
      progress = sprintf('\r%03.1f%%', progress)
      cat(progress)
      sourced = F
      lastSource = NA ## index, not uid
      lastInit = NA ## index, not uid
      reqState = factor('inactive', levels=reqStates)
      for(reqUID in uidsByReq[[req]]){
        reqIndex = x[J(reqUID), which=T]
        if(!reqIndex){
          stop("missing UID: ", reqUID, " in x!\n")
        }
        if(reqState == 'inactive'){
          if(x[reqIndex, name] %in% MPI_Req_inits){
            lastInit = reqIndex
            reqState = factor('init', levels=reqStates)
          } else if(x[reqIndex, name] %in% MPI_Req_sources){
            lastSource = reqIndex
            reqState = factor('active', levels=reqStates)
          } else {
            reqWarning()
            break
          }
        } else if(reqState == 'init'){
          if(x[reqIndex, name] %in% MPI_Req_starts){
            lastSource = reqIndex
            reqState = factor('active', levels=reqStates)
            ##! MPI_Startall has multiple reqs
            brefs[[reqIndex]] = union(brefs[[reqIndex]], lastInit)
            frefs[[lastInit]] = union(frefs[[lastInit]], reqIndex)
          } else if(x[reqIndex, name] %in% MPI_Req_frees){
            reqState = factor('inactive')
          } else {
            reqWarning()
            break
          }
        } else if(reqState == 'active'){
          if(x[reqIndex, name] %in% MPI_Req_sinks){
### tests are conditional sinks; unsuccessful tests don't record request IDs
            ## wait or successful test
            reqState = factor('null', levels=reqStates)
            brefs[[reqIndex]] = union(brefs[[reqIndex]], lastSource)
            frefs[[lastSource]] = union(frefs[[lastSource]], reqIndex)
          } else if(x[reqIndex, name] %in% MPI_Req_frees){
            reqState = factor('inactive', levels=reqStates)
          } else {
            reqWarning()
            break
          }
        } else if(reqState == 'null'){
          if(x[reqIndex, name] %in% c(MPI_Req_sources, MPI_Req_starts)){
            lastSource = reqIndex
            reqState = factor('active', levels=reqStates)
            brefs[[reqIndex]] = union(brefs[[reqIndex]], lastInit)
            if(x[reqIndex, name] %in% MPI_Req_starts)
              frefs[[lastInit]] = union(frefs[[lastInit]], reqIndex)
          } else if(x[reqIndex, name] %in% MPI_Req_frees){
            reqState = factor('inactive')
            ##!@todo reqIndex is 0 for nekbone; a uid is missing, so line 370 assigns 0 to 
          } else if(x[reqIndex, name] %in% MPI_Req_inits){
            lastInit = reqIndex
            reqState = factor('init')
          } else {
            reqWarning()
            break
          }
        }
        offset = offset + 1
      }
    }
    cat('\r100.0%\n')
    ##! these are indices, but need to be uids
    sel = !sapply(brefs, is.null)
    ## x key must be uid
    brefs[sel] = lapply(brefs[sel], function(b) x[b, uid])
    brefs[!sel] = as.numeric(NA)
    x$bref = brefs

    sel = !sapply(frefs, is.null)
    frefs[sel] = lapply(frefs[sel], function(b) x[b, uid])
    frefs[!sel] = as.numeric(NA)
    x$fref = frefs
    rm(brefs, frefs)
  } ## if(length(requests))

  if(!debug)
    x[,reqs:=NULL]

  ## remove src and tag from req sinks, convert back from list
  x[name %in% MPI_Req_sinks, src := as.integer(NA)]
  x[, src := unlist(src)]

  x[name %in% MPI_Req_sinks, tag := as.integer(NA)]
  x[, tag := unlist(tag)]

  if(writeReplay){
    #!@todo make this a function
    cat('Rank', rank, 'writing replay\n')
    ##!@todo the replay fscanf is going to trip on NAs
    cols = names(colClasses)
    cols = grep('_w$', cols, v=T, invert=T)
    y = data.table::copy(x[, cols, with=F][order(start)])
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

preDeps = function(x, uidsByReq, ...){
  ranks = as.integer(names(x))
  
  x = mcmapply(.deps, x, ranks, uidsByReq, ..., SIMPLIFY=F)
  return(x)
}

deps = function(x){
  if(debug)
    cat('Merging tables\n')
  ## merge tables
  commTable = .rbindlist(lapply(x, '[[', 'comms'))
 
  x = .rbindlist(lapply(x, '[[', 'runtimes'))[order(start)]
  ranks = sort(unique(x$rank))
  
  if(debug)
    cat('Unifying communicators\n')
  if(nrow(commTable) > 0){
    ## unify communicators
    commTable$done = F
    cid = 1
    
    commList = list(new_MPI_COMM_WORLD) ## -1 is new MPI_COMM_WORLD
    while(length(commList)){
      comm = commList[[1]]
      commList = tail(commList, -1)
      setkey(commTable, rank, source, done)
      sel = commTable[parentComm == comm & !done, which=T]

      d = commTable[sel,list(source=min(source)),by=rank]
      d$done = F
      setkey(d)
      sel = commTable[d, which=T]
      sel = intersect(sel, commTable[childComm != new_MPI_COMM_NULL, which=T])

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

      ##! this is applied on a per-rank basis
      f = function(row)
        x[uid >= row$source &
          uid <= row$sink &
          comm == row$childComm &
          rank == row$rank,
          comm := row$unifiedComm]
      rowApply(translate, f)
      rm(translate)
    
      ## if more entries in current comm, continue
      if(nrow(commTable[parentComm == comm & !done]))
        commList = c(comm, commList)
      
      ## get next comm (must be a child comm by definition)
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

  x[name %in% MPI_Inits, 'vertex'] = 1
  x[name == 'MPI_Finalize', 'vertex'] = 2

  if(debug)
    cat('Collectives\n')
  collectives = intersect(unique(x$name), setdiff(MPI_collectives, c('MPI_Init','MPI_Finalize')))

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
  collectives = 1:(vid-1)

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
  return(list(runtimes = x, comms = commTable, collectives=collectives))
}

messageDeps = function(x){
  ## match inter-rank messages. Messages (excluding MPI_ANY_SOURCE and
  ## MPI_ANY_TAG) are uniquely identified by the following: source,
  ## destination, size, tag, communicator, and order.  For
  ## MPI_ANY_SOURCE and MPI_ANY_TAG, we record the tag and source at
  ## run time.

  commTable = x$comms
  x = x$runtimes
  messageCols = c('src', 'dest', 'size', 'tag', 'comm')
  nonMessageCols = setdiff(names(x), messageCols)

  mids =
    unique(x[name %in% c(MPI_Sends, MPI_Recvs, MPI_Req_inits),
             messageCols, with=F])
  sel = complete.cases(mids)
  if(length(which(!sel))){
    cat('Removing', length(which(!sel)), 'message IDs:\n')
    print(mids[!sel])
    mids = mids[sel]
  }

  if(nrow(mids) < 1){
    cat('no messages\n')
    return(list(runtimes=x, messages=NULL))
  }

  xStarts = x[name %in% MPI_Req_starts]
  setkey(xStarts, uid)
  setkey(x, src, dest, size, tag, comm)
  canceledUIDs = x[name == 'MPI_Cancel', uid]
  
  f_noSideEffects = function(s, r){
    ## handle Irecvs by forwarding to corresponding MPI_Wait()s
    if(!is.na(r$fref)){ ## receive is waited for later
      dest = unlist(r$fref) ## uid
      o_dest = r$uid ## uid
      dest_vertex = as.numeric(NA)
      o_dest_vertex = r$vertex
    } else {
      dest = r$uid ## uid
      o_dest = dest ## uid
      dest_vertex = r$vertex
      o_dest_vertex = as.numeric(NA)
    }
    src = s$uid ## uid
    o_src = src ## uid
    tryCatch(
      result <-
      data.frame(o_src=o_src, src=src, o_dest=o_dest, dest=dest,
                 size=s$size,
                 src_vertex=s$vertex, dest_vertex=dest_vertex,
                 o_dest_vertex=o_dest_vertex,
                 src_rank=s$rank, dest_rank=r$rank, tag=s$tag),
      error = function(e) {
        print(s)
        print(r)
        stop('error in f_noSideEffects: ')
      }
      )
    return(result)
  }

  f = function(mid){
    matching = x[mids[mid]]

    ## ignore canceled requests
    matching = matching[is.na(fref) | !fref %in% canceledUIDs]
    ## ignore requests not waited for
    matching = matching[!name %in% c(MPI_Req_sources, MPI_Req_inits) | !is.na(fref)]
    if(nrow(matching) < 1){
      cat('Message ID', mid, 'of', nrow(mids), ':', 'no messages\n')
      return(data.frame())
    }

    setkey(matching, uid)
    sends = matching[name %in% MPI_Sends]
    recvs = matching[name %in% MPI_Recvs]
    
    if(nrow(xStarts)){
      sendInit_frefs = unlist(matching[name %in% MPI_Send_inits, fref])
      if(length(sendInit_frefs)){
        sendStarts = xStarts[J(sendInit_frefs), nonMessageCols, with=F]
###!@ todo all messageCols should be identical, so something like this
###!should work: .rbindlist(rep.int(mids[mid], nrow(sendStarts))))
        sendStarts =
          cbind(sendStarts, 
                matching[J(unlist(sendStarts[, bref])), messageCols, with=F])
        setcolorder(sendStarts, names(x))
        sends = rbind(sends, sendStarts)
        rm(sendStarts, sendInit_frefs)
      }
      recvInit_frefs = unlist(matching[name %in% MPI_Recv_inits, fref])
      if(length(recvInit_frefs)){
        recvStarts = xStarts[J(recvInit_frefs), nonMessageCols, with=F]
        recvStarts =
          cbind(recvStarts, matching[J(unlist(recvStarts[, bref])), messageCols, with=F])
        setcolorder(recvStarts, names(x))
        recvs = rbind(recvs, recvStarts)
        rm(recvStarts, recvInit_frefs)
      }
    }
    sends = sends[order(uid)]
    recvs = recvs[order(uid)]
    
    if(debug && !(mid %% 100))
      cat('Message ID', mid, 'of', nrow(mids), ':', nrow(sends),
          'messages\n')
    ## if a receive has a request, the send dependency leads to the
    ## matching wait, test, or free
    if(nrow(sends) != nrow(recvs)){
      errMsg = paste('mismatched send-receive:\n',
        paste(names(mids), collapse='\t'), '\n',
        paste(mids[mid], collapse='\t'))
      stop(errMsg)
    }

    result =
      lapply(1:nrow(sends), function(row)
             f_noSideEffects(sends[row], recvs[row]))
    result = .rbindlist(result)
    return(result)
  }

  ## srDeps holds source and destination UIDs for matching messages
  srDeps = mclapply(1:nrow(mids), f)
  srDeps = srDeps[sapply(srDeps, nrow) > 0]
  ## this is slower than rbindlist, but doesn't segfault. rbindlist
  ## makes shitty data tables?
  srDeps = .rbindlist(srDeps)
  srDeps = as.data.table(lapply(srDeps, unlist))

  ## find dest vertices for nonblocking receives
  setkey(x, uid)
  missingVertices = srDeps[is.na(dest_vertex), which=T]
  if(length(missingVertices)){
    dest_uids = srDeps[missingVertices, dest]
    srDeps[missingVertices, dest_vertex := x[J(dest_uids), vertex]]
  }
  
  if(debug)
    cat('Done finding message matches\n')


  setkey(x, uid)
  ## find indices with multiple references
  setkey(srDeps, dest)
  destTable = table(srDeps[, dest])
  multDests = as.numeric(names(destTable[destTable > 1]))

  ## complete the indices with single references
  singleDests = as.numeric(names(destTable[destTable == 1]))
  deps = vector('list', nrow(x))
  if(length(singleDests))
    deps[x[J(singleDests), which=T]] = srDeps[J(singleDests)]$src ## uid

  ## complete the indices with multiple references
  if(length(multDests)){
    sel = x[J(multDests), which=T]
    deps[sel] =
      mapply(c, deps[sel],
             lapply(multDests,
                    function(d) srDeps[J(d)]$src), SIMPLIFY=F)
  }
  x$deps = deps
  rm(deps)
    
  if(debug)
    cat('Done matching messages\n')

  ## "messages" columns are UIDs
  return(list(runtimes=x, messages=srDeps))
}

### vertices and edges with attributes.
tableToGraph = function(x, assignments, messages, saveGraph=T, path='.'){

  #host = unique(sub('[[:digit:]]+', '', assignments$hostname))

  if(debug)
    cat('Computation edges\n')
  startTime = Sys.time()
  ## replace computation vertices with weighted edges
  setkey(x, rank, uid)
  f = function(r){
    sel = x[rank == r & is.na(name), which=T]
    result =
      cbind(x[sel, list(weight=duration,
                        power=pkg_w+dram_w,s_uid=uid,flags)],
            x[sel-1,list(src=vertex)],
            x[sel+1,list(dest=vertex, d_uid=uid)])
    result[, rank:=r]
  }
  compEdges = .rbindlist(mclapply(unique(x[, rank]), f))
  
  ## delete computation rows
  x = x[!is.na(name)]

  cat('Comp edges time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()

  ##!@todo decide how to handle collectives (decompose, single vertex, etc.)

  ## For the vertex frame, column 1 is the vertex name.
  ##vertices = x[, list(name, size, dest, src, tag, comm, hash, vertex)]
  vertices = x[, list(name, rank, hash, vertex, reqs)]
  vertices$name = paste(vertices$name, vertices$rank)
  setnames(vertices, 'name', 'label')
  vertices$rank = NULL
  setcolorder(vertices, c('vertex', setdiff(names(vertices), 'vertex')))
  setkey(vertices, vertex)
  vertices = unique(vertices)

  ## For the edge list frame, column 1 is the source and column 2 is the dest
  setkey(x, uid)
  setkey(assignments, rank)
  ## assume all ranks run on a cluster with a single hostname prefix
  uhost = unique(sub('[[:digit:]]+', '', assignments$hostname))
  
### compute message source and dest vertices with weight (latency)

### o_dest and dest uids should be on the same rank anyway
  
  if(debug)
    cat('Communication edges\n')
  ## add communication edges
  edges = compEdges[,list(src,dest,weight,s_uid,d_uid)]
  if(!is.null(messages) && nrow(messages)){
    setkey(messages, src_rank)
    messages$src_host = assignments[messages, hostname]

    setkey(messages, dest_rank)
    messages$dest_host = assignments[messages, hostname]

    messages[src_host == dest_host,
             weight := selfLatency[[uhost]](size), by=list(o_src, o_dest)]
    messages[src_host != dest_host,
             weight := latency[[uhost]](size), by=list(o_src, o_dest)]

    messageEdges =
      messages[, list(size, src=src_vertex, dest=dest_vertex,
                      o_dest=o_dest_vertex,
                      rank=dest_rank, s_rank=src_rank,
                      s_uid=o_src, ##!@todo is this correct for MPI_Send_init/MPI_Start?
                      d_uid=dest,
                      o_d_uid=o_dest,
                      weight, tag)]
    rm(messages)
    commonNames = intersect(names(edges), names(messageEdges))
    messageOnlyNames = setdiff(names(messageEdges), commonNames)
    for(name in messageOnlyNames)
      edges[, name := as.numeric(NA), with=F]
    edges = rbind(edges, messageEdges, use.names=T)
  } else
    messageEdges = NA
    
  cat('Comm edges time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  if(debug)
    cat('Graph object\n')
  g = graph.data.frame(edges, vertices=vertices)
  rm(edges)
  if(saveGraph){
    graphFile = file.path(path, 'graph.dot')
    write.graph(g, file=graphFile, format='dot')
    system(paste('gzip -f ', graphFile, sep=''), wait=F)
  }
  cat('Graph object time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  return(list(graph=g, messageEdges = messageEdges, compEdges = compEdges,
              vertices=vertices))
}

tableToMarkov = function(x, rank=0, path='.'){
  .rank = rank
  setkey(x, rank)
  x = x[rank == .rank]
  setkey(x, hash)
  y = data.table::copy(x)
  setkey(y, uid)
  hashes = unique(x[!is.na(name) & name != 'MPI_Finalize']$hash)
  vertices =  x[!is.na(name), list(label=unique(name)), by=hash]
  edges =
    .rbindlist(
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

## Populate power for each task. For now, overwrite actual
## measurements.
##
## Power numbers depend on how many sockets are used by each rank.
modelPower = function(x, assignments){
  host = unique(sub('[[:digit:]]+', '', assignments$hostname))
  cpu = cpuMap[[host]]
  fActivePower = activePower[[cpu]]
  idlePower = idlePower[[cpu]]

  x[, pkg_w:=NULL]
  x[, pp0_w:=NULL]
  x[, dram_w:=NULL]

  hostCount = assignments[, list(count=nrow(.SD)), by=list(hostname)]

  sharedHosts = hostCount[count > 1]$hostname
  if(length(sharedHosts)){
    cat('Power estimation with host sharing is not yet supported\n')
    x$total_power = as.numeric(NA)
    return(x)
    
    ## check for socket sharing
    socketCount =
      assignments[, list(count=nrow(.SD)), by=list(hostname, socket)]
    if(any(socketCount$count > 1)){
      cat('Cannot determine power; multiple ranks share a socket.')
      x$total_power = as.numeric(NA)
      return(x)
    }

    ##!@todo omp_sockets is in string format; convert to integer list for comparison
    omp_socketSharing = sapply(sharedHosts, function(h){
      sharedSockets =
        do.call(intersect, assignments[hostname == h]$omp_sockets)
      return(length(sharedSockets) > 0)
    })
    if(any(omp_socketSharing)){
      cat('Cannot determine power; multiple ranks share a socket.')
      x$total_power = as.numeric(NA)
      return(x)
    }
  }
  ## determine number of sockets in use by each rank

  ##!@todo get thread count, serial time, and parallel time from
  ##instrumentation.
  
  ##x[, total_power := fActivePower(threads, cpu_freq)]
  
  x
}

shortStats = function(x, thresh=.001){
  shortTaskRatio = nrow(x[duration < thresh])/nrow(x)
  shortTimeRatio = sum(x[duration < thresh, duration])/sum(x$duration)
  cat(shortTaskRatio * 100, '% of tasks <', thresh, 's\n')
  cat(shortTimeRatio * 100, '% of time in short tasks\n')
}

run = function(path='.', saveResult=F, name='merged.Rsave', noReturn=F){
  startTime = Sys.time()
  cat('Reading ', path,' (', getwd(), ')\n')
  a = readAll(path)
  cat('Read time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  assignments = a$assignments
  uidsByReq = a$uidsByReq
  newComms = a$newComms
  b = preDeps(a$runtimes, uidsByReq, newComms, path=path)
  cat('preDeps time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  rm(a)
  b = deps(b)
  cat('deps time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  comms = b$comms
  collectives = b$collectives ## vertex IDs
  b2 = messageDeps(b)
  cat('messageDeps time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  messages = b2$messages
  b2 = b2$runtimes
  rm(b)
  ##if(!saveResult)
  g = tableToGraph(b2, assignments=assignments, messages=messages, path=path)
  ##else
  ##  g = NA
  cat('tableToGraph time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  startTime = Sys.time()
  #tableToMarkov(b2, path=path)
  #cat('Markov time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  #startTime = Sys.time()
  
  result =
    list(##runtimes = b2,
         messageEdges = g$messageEdges,
         compEdges = g$compEdges,

         ##!@todo some vertices are missing...
         vertices = g$vertices,
         
         assignments = assignments,
         comms = comms,
         globals=globals,
         collectives = collectives)
  if(saveResult){
    with(result, save(list=ls(), file=file.path(path,name)))
    cat('save time: ', difftime(Sys.time(), startTime, units='secs'), 's\n')
  } else
    result$graph = g$graph
  if(noReturn)
    return(NULL)
  else {
      result$runtimes = b2
      return(result)
  }
}
