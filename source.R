# collapse.rows()
collapse.rows <- function(expr, probe.col, gene.col, data.table=F, method=c("maxMean", "minMean", "colMean", "colMedian")){
  if(length(grep('data.table', installed.packages())) == 0){
    install.packages('data.table')
    require(data.table)
  }else if(length(grep('data.table', search())) == 0){
    suppressPackageStartupMessages(require(data.table))
  }
  
  if (probe.col == "rownames"){
    expr <- data.table(expr, keep.rownames=T)
    setnames(expr, "rn", "rownames")
  }else{
    expr <- data.table(expr)
  }
  
  if(method=="maxMean" | method=="minMean"){ 
    expr[, rowmean := rowMeans(.SD[, !c(probe.col, gene.col), with=F])]
    if(method=="maxMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[1], by=gene.col][, rowmean:=NULL]
    }
    else if(method=="minMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[.N], by=gene.col][, rowmean:=NULL]
    }
  }
  else if(method=="colMean"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], mean), by=gene.col]
  }
  else if(method=="colMedian"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], median), by=gene.col]   
  }
  else stop("method must be 'maxMean', 'minMean', 'colMean' or 'colMedian'\n")
  
  if(!data.table){
    return(data.frame(res))
  }else{ return(res[]) }
}