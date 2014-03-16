#!/usr/bin/env Rscript

library('data.table')
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
  result = rbindlist(result)
  return(result)
}

deps = function(x){
  x$vertex = as.numeric(NA)

  ## init and finalize
  init = which(x$name == 'MPI_Init' | x$name == 'MPI_Init_thread')
  finalize = which(x$name == 'MPI_Finalize')
  vid = 3

  ## g = data.table(name='MPI_Init', succ=as.integer(NA), vertex=1)
  x[init, 'vertex'] = 1
  
  ## g = rbindlist(list(g, list(name='MPI_Finalize', succ=NA, vertex=2)))
  x[finalize, 'vertex'] = 2

  ##!@todo parse list of communicators and their members
  ## sequential dependencies per rank
  ranks = unique(x$rank)
  for(rank in ranks){
    rows = which(x$rank == rank)
    x$deps[tail(rows, -1)] = head(rows, -1)
  }

  ## collectives
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

  ## handle blocking receives and waits, and nonblocking successful tests
  
  return(list(table=x, graph=g))
}

