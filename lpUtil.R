require('rjson')
require('data.table')
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
        b
      })
    } else
      arrayVars = NULL
    if(length(singleSel)){
      singleVars = lapply(b[singleSel], as.data.table)
    } else
      singleVars = NULL
    c(arrayVars, singleVars)
  }
  a$Solution[[2]]$Variable = f(a$Solution[[2]]$Variable)
  a$Solution[[2]]$Constraint = f(a$Solution[[2]]$Constraint)
  
  a
}

writeSolve = function(){
  confName = writeSlices(reduced)
  
}
