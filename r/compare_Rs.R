
#library(socialmixr)







compare_Rs = function(eg_props_s, breaks = c(0,5,12,18,30,40,50,60,70,Inf), tranvec = NULL, suscvec = NULL)
  {

  if(is.null(tranvec)){ tranvec = rep(1, length(breaks)-1)}
  if(is.null(suscvec)){ suscvec = rep(1, length(breaks)-1)}
  
  
  q = outer(tranvec, suscvec)
  
  
  get_R = function(eg, breaks = breaks, q = q){
    R = max(eigen(matrix(eg$aug_mean_sym, nrow = length(breaks)-1)*q)$values)
    return(R)
    
  }
  
  
  Rs = sapply(eg_props_s, FUN = get_R, breaks = breaks, q = q)
  Rs
  
  # breaks[length(breaks)] <- 120
  # data(polymod)
  # cm <- contact_matrix(polymod, countries = "United Kingdom", age.limits = breaks, symmetric = TRUE)
  # 
  # ngm_polymod = cm$matrix * q
  # 
  # pmR <- max(eigen(ngm_polymod)$values)
  # 
  # Rs/pmR
  # 
  # 
  # rvals <- rnorm(1000, 2.6, sd = 0.54)
  # 
  # Rvals = lapply(Rs/pmR, FUN = function(X){X * rvals})
  # Rvals

}


compare_Rs_stan= function( sym_mats_s, breaks = c(0,5,12,18,30,40,50,60,70,Inf), weeks_range = 1:33, tranvec = NULL, suscvec = NULL){
  
  if(is.null(tranvec)){ tranvec = rep(1, length(breaks)-1)}
  if(is.null(suscvec)){ suscvec = rep(1, length(breaks)-1)}
  
  
  q = outer(tranvec, suscvec)
  
  levs <- unique(unlist(cut(seq(0,120),breaks, right=FALSE), use.names = FALSE))
  eg = expand.grid(sort(levs),sort(levs))

  get_Rs = function(n, sym_mats_s_, eg_=data.table(eg),  breaks_, q_){
    eigs = c()  
    samps = length(sym_mats_s_[[n]])/(length(breaks)-1)^2
    for(i in 1:samps){
      means = sym_mats_s_[[n]][,i]
      eg_[,aug_mean_sym:=means]
      eigs  = append(eigs,max(Re(eigen(matrix(eg_$aug_mean_sym, nrow = length(breaks_)-1)*q_)$values)))
    }
  return(eigs)
  }
  
  eigss = lapply(weeks_range, get_Rs, breaks_ = breaks, q_ = q, sym_mats_s_=sym_mats_s)
  
  matrix(unlist(eigss), nrow = length(eigss), byrow = T)
}  


