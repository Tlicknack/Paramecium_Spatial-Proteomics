NAprocessPCP_SUP = function(msn){
  msn_sup = msn[which(msn[,ncol(msn)] > 0),]
  for(i in ((ncol(msn)-11):(ncol(msn)-1)) ){
    msn_sup = msn_sup[which(is.na(msn_sup[,i])),]
  }
  return(msn_sup)
}