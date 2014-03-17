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
  a[size == -1,]$size = NA
  a[dest == -1,]$dest = NA
  a[src == -1,]$src = NA
  a[is.na(size)][tag == -1]$tag = NA
  a$reqs = strsplit(a$reqs,',')
  a[comm == '0x0',comm:='']
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
.deps = function(x){
  x$vertex = as.numeric(NA)
  x$deps = vector('list', nrow(x))
  x$brefs = vector('list', nrow(x))
  x$fref = as.numeric(NA)

  ##!@todo parse list of communicators and their members
  ## sequential dependencies per rank
  ranks = unique(x$rank)
  if(length(ranks) > 1)
    stop('Too many ranks in deps()!')

  rows = 1:nrow(x)
  x$deps[tail(rows, -1)] = head(rows, -1)

  ## get backrefs for blocking receives and waits,
  ## nonblocking successful tests, and request frees
  if(debug)
    cat('Rank', ranks, 'requests\n')
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
        cat('Rank', ranks, 'sink request', sinkIndex, 'of', length(sinkReqs), '\n')
      reqFound =
        reqs[sapply(x[reqs]$reqs, function(x) req %in% x)]
###      cat('Rank', ranks, 'reqFound:', reqFound, '\n')
      sourced = F
      for(i in reqFound){
        if(x[i]$name %in% MPI_req_sources){
          sourced = T
          last = i
        } else if(sourced){ ## sink
###          if(debug)
###            cat('Rank', ranks, 'link:', last, 'to', i, '\n')
          ##!@todo fix multiple assignment
          x[i]$brefs = union(unlist(x[i]$brefs), last)
          x[last]$fref = i
          sourced = F
        }
      }
      sinkIndex = sinkIndex + 1
    }
  }

  return(x)
}

deps = function(x){
  ##if(debug)
  ##  x = lapply(x, .deps)
  ##else
    x = mclapply(x, .deps)

  ## merge tables
  ##!@todo adjust indices for dpes, bref, and fref 
  x = rbindlist(x)[order(start)]
  
  ranks = unique(x$rank)
  
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
  
  return(x)
  ##return(list(table=x, graph=g))
}
