#!/usr/bin/env Rscript

setwd(commandArgs(trailingOnly=T))

#!@todo check if result file exists; if so, quit

source('../read.R')
a = run()
s = list()
wd = tail(strsplit(getwd(), '/')[[1]], 1)
collectives = setdiff(MPI_collectives, c('MPI_Init', 'MPI_Finalize'))
##s$offenders = a$runtimes[rank==0 & name %in% collectives & duration > .1]
s$offenders = a$runtimes[name %in% collectives, {s = .SD[, list(maxDur=max(duration), arrSpread = diff(range(start)))]; list(offending= s$maxDur >= 2* s$arrSpread && s$maxDur > .1)}, by=vertex]
s$date = wd
s$ranks = unique(a$runtimes$ranks)

if(any(s$offenders$offending))
  save(s, file=paste('../', wd, '.long_collectives.dat', sep=''))
