NAprocessPCP_MAC = function(msn){
  msn_mac = msn[which(msn[,ncol(msn)-11] > 0),]
  for(i in (ncol(msn)-10):(ncol(msn)) ){
    msn_mac = msn_mac[which(is.na(msn_mac[,i])),]
  }
  return(msn_mac)
}
