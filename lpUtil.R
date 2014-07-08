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

reconcileLP = function(resultFile, timesliceFile){
  result = readLP(resultFile)
  load(timesliceFile)
  list(result=result, slice=slice)
}

timeStr = '[0-9]+[.][0-9]+'
readCommandResults = function(command){
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
    times = sub('p.*w[.]results$', '',
      sub(paste(command, '_', sep=''), '', resultFiles))
    result = mcmapply(reconcileLP, resultFiles, timesliceFiles, SIMPLIFY=F)
    names(result) = times
    result
  }
  nnapply(powerLimits, f)
}

files = list.files(pattern='.*[.]results$')
commands = unique(sub(paste('_', timeStr, '[.]p.*w[.]results', sep=''),'',files))
results = nnapply(commands, readCommandResults)
