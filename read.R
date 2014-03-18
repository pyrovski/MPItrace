#!/usr/bin/env Rscript

debug=T

library('data.table')
library('parallel')
source('~/local/bin/pbutils.R')

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
  'MPI_File_write_at_all')

MPI_req_sinks =
  c('MPI_Wait', 'MPI_Waitall', 'MPI_Waitsome',
    'MPI_Test', 'MPI_Testall', 'MPI_Testsome',
    'MPI_Request_free')

MPI_req_sources =
  c('MPI_Isend','MPI_Irecv','MPI_Irsend','MPI_Bsend_init',
    'MPI_Rsend_init','MPI_Send_init','MPI_Ssend_init',
    'MPI_Recv_init')

readGlobal = function(filename = "glog.dat"){
  eval(parse(filename), envir=.GlobalEnv)
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
  files = sort(list.files(path,'runtime.*.dat'))
  ranks = as.numeric(t(unname(as.data.frame(strsplit(files,'[.]'))))[,2])
  result = lapply(files, readRuntime)
  names(result) = ranks
  result = napply(result, function(x, name){x$rank = as.numeric(name);x})
  ## result = rbindlist(result)
  return(result)
}

## single rank only!
.deps = function(x, maxRank){
  ranks = unique(x$rank)
  if(length(ranks) > 1)
    stop('Too many ranks in deps()!')

  rank = ranks
  rm(ranks)

  x$vertex = as.numeric(NA)

  # numeric for now, list later
  x$deps = as.numeric(NA)
  x$brefs = vector('list', nrow(x))
  x$fref = as.numeric(NA)
  x$uid = (1:nrow(x))+(rank+1)/(maxRank+1)

  ##!@todo parse list of communicators and their members

  ## sequential dependencies

  rows = 1:nrow(x)
  x$deps[tail(rows, -1)] = x$uid[head(rows, -1)]

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
          ## if(debug){
          ##   cat('Rank', rank, 'link:', last, 'to', i, '\n')
          ##   cat('Rank', rank, 'brefs before adding', x$uid[last],':', unlist(x[i]$brefs), '\n')
          ## }
          ##!@todo fix multiple assignment
          x$brefs[[i]] = union(unlist(x[i]$brefs), x$uid[last])
          x[last]$fref = x$uid[i]
          ## if(debug){
          ##   cat('Rank', rank, 'brefs after:', unlist(x[i]$brefs), '\n')
          ## }
          sourced = F
        }
      }
      sinkIndex = sinkIndex + 1
    }
  }
  return(x)
}

deps = function(x){
  ranks = sapply(x, function(x) unique(x$rank))

  ## if(debug)
  ##   x = lapply(x, .deps, max(ranks))
  ## else
    x = mclapply(x, .deps, maxRank=max(ranks))

  ## merge tables
  ##!@todo order by topological sort after all dependecies done
  x = rbindlist(x)[order(start)]
  
  ## remap uids
  newUIDs = as.list(1:nrow(x))
  names(newUIDs) = x$uid
  x$uid = unlist(newUIDs[as.character(x$uid)])

  sel = x[!is.na(fref), which=T]
  x$fref[sel] = unlist(newUIDs[as.character(x$fref[sel])])

  sel = x[!is.na(deps), which=T]
  x$deps[sel] = newUIDs[as.character(x$deps[sel])]
  
  sel = which(!sapply(x$brefs, is.null))
  ##!@todo speed this up
  x$brefs[sel] = sapply(x$brefs[sel], function(x)unlist(unname(newUIDs[sapply(x,as.character)])))

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
  
  setkey(x, name, comm)
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

  ##!@todo match inter-rank messages

  ##!@todo add vertices for intra-rank stuff
  ##!@todo decide how to handle collectives (decompose, etc.)
  
  return(x)
  ##return(list(table=x, graph=g))
}
