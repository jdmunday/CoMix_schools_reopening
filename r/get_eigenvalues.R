


get_eigenvalues = function(fit_type, samples, cap, week_range,  nwks=2, breaks=c(0,5,12,18,30,40,50,60,70,Inf), cms_all = NULL, suscvec=NULL, tranvec=NULL)
  {

  if (is.null(cms_all)){
    
    # initialise empty containers --------------------------------------------------
    
    cms_all = list()
    cms = list()
    
    # Load in pre-calculated contact matrices--------------------------------------------------
    
    for(week in week_range){
      cms[[week]] = qs::qread(paste0('data/contact_matrices/', fit_type, samples,  '_ngrps', length(breaks) - 1,'_cap', cap, '_nwks', nwks, '_sr', week,'_scms.qs'))
      }
  
  
    cms_all = cms
  }
  else{
    print('using provided cms')
  }

  
  # Get susceptibility vector --------------------------------------------------
  
  
  
  # Calculate Rs from eigs relative to polymod --------------------------------------------------
  
  eigs <- compare_Rs_stan(cms_all, breaks, weeks_range = week_range, suscvec = suscvec, tranvec = tranvec)
  return(eigs)
}

calculate_eigs_and_rats = function(eig_stats_ = eig_stats, eig_rats_ = eig_rats, eig_frame_= eig_frame, suscvec_ = rep(1,9), tranvec_=rep(1,9), suscname='1. Equal', all_cms_NP, all_cms_NS){
  
  
  eigs_LK2NS = get_eigenvalues('bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NS, suscvec = suscvec_, tranvec = tranvec_)
  eigs_LK2NP = get_eigenvalues('bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NP, suscvec = suscvec_, tranvec = tranvec_)
  eigs_lK3AS = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_, tranvec = tranvec_)
  eigs_lK2AS = get_eigenvalues('bs', 1000, 50, 33:33, nwks=5,  suscvec = suscvec_, tranvec = tranvec_)
  
  
  periods = c('Primary', 'Secondary', 'Both', 'None')
  suscrulelist = rep(suscname, 4)
  i = 1
  for(rats in list(eigs_LK2NP, eigs_LK2NP, eigs_lK2AS, eigs_lK3AS)){
    eig_stats_ = rbind(eig_stats_, data.table(
      
      compare = periods[i], 
      suscrule = suscrulelist[i],
      compare_num = i,
      rat05 = quantile(ecdf(rats), 0.05),
      rat50 = quantile(ecdf(rats), 0.5),
      rat95 = quantile(ecdf(rats), 0.95),
      rat25 = quantile(ecdf(rats), 0.25),
      rat75 = quantile(ecdf(rats), 0.75),
      rat10 = quantile(ecdf(rats), 0.10),
      rat90 = quantile(ecdf(rats), 0.9)
      
      
      
    ))
    i = i + 1
  }
  
  l2NPRAT = eigs_LK2NP / eigs_lK3AS
  l2NSRAT = eigs_LK2NS / eigs_lK3AS
  l2ASRAT = eigs_lK2AS / eigs_lK3AS
  
  
  periods = c('Primary', 'Secondary', 'Both')
  suscrulelist = rep(suscname, 3)
  i = 1
  for(rats in list(l2NSRAT, l2NPRAT, l2ASRAT)){
    eig_rats_ = rbind(eig_rats_, data.table(
      
      compare = periods[i], 
      suscrule = suscrulelist[i],
      compare_num = i,
      rat05 = quantile(ecdf(rats), 0.05),
      rat50 = quantile(ecdf(rats), 0.5),
      rat95 = quantile(ecdf(rats), 0.95),
      rat25 = quantile(ecdf(rats), 0.25),
      rat75 = quantile(ecdf(rats), 0.75),
      rat10 = quantile(ecdf(rats), 0.10),
      rat90 = quantile(ecdf(rats), 0.9)
      
      
      
    ))
    i = i + 1
  }
  
  i = 1
  periods = c('Primary', 'Secondary', 'Both')
  suscrulelist = rep(suscname, 3)
  for(rats in list(l2NSRAT, l2NPRAT, l2ASRAT)){
    eig_frame_ = rbind(eig_frame_, data.table(
      
      Period = rep(periods[i], length(rats)), 
      Ratio = as.vector(rats),
      suscrule = rep(suscname, length(rats))
      
    ))
    i = i + 1
  }
  
  list(eig_stats_, eig_rats_, eig_frame_)
  
}



