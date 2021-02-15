


get_all_cms_vals = function(area_name, fit_type, samples, cap, week_range,  nwks=2, breaks=c(0,5,12,18,30,40,50,60,70,Inf), cms_all = NULL, suscvec=NULL, tranvec=NULL, usefirst=1000)
{
  
  
  if (is.null(cms_all)){
    
    # initialise empty containers --------------------------------------------------
    
    cms_all = list()
    cms = list()
    
    # Load in pre-calculated contact matrices--------------------------------------------------
    
    for(week in week_range){
      cms[[week]] = qs::qread(paste0('data/contact_matrices/', fit_type, samples,  '_ngrps', length(breaks) - 1,'_cap', cap, '_nwks', nwks, '_sr', week,'_scms.qs'))[,1:usefirst]
    }
    
    
    cms_all = cms
  }
  else{
    print('using provided cms')
  }
  
  
  # Get susceptibility vector --------------------------------------------------
  
  # Calculate Rs from eigs relative to polymod --------------------------------------------------
  
  cms_all
  
}


get_all_egs = function(
  filenames = list( "data/contact_matrices/bs1000_ngrps9_cap50_nwks5_sr33_scms.qs", 
                    "data/contact_matrices/bs1000_ngrps9_cap50_nwks5_sr41_scms.qs"),
  periods = c('Lockdown 2', 'Lockdown 3')
){
  
  all_egs = data.table()
  i = 1
  for( fn in filenames){
    print(fn)
    

    
    #dts = parts[!panel %in% c("C", "D") & survey_round %in% weeks, .(start_date = min(date), end_date = max(date))]
    
    cms = qs::qread(fn)
    levs <- unique(unlist(cut(seq(0,120),breaks, right=FALSE), use.names = FALSE))
    age_lab <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")
    
    # Get columns of age-groups to put into mapply
    eg = expand.grid(sort(levs),sort(levs))
    
    
    
    eg = data.table(eg)
    eg[, Var1 := factor(Var1, levels = levs, labels = age_lab)]
    eg[, Var2 := factor(Var2, levels = levs, labels = age_lab)]
    
    eg[,cms := rowMeans(cms)]
    
    eg[,sr:=periods[i]]
    
    #eg[, daterange := paste0(dts$start_date, ' - ', dts$end_date)]
    
    
    all_egs = rbind(all_egs, eg)
    
    i = i + 1
  }
  
  all_egs
}



augment_cms = function(
  filenames = list( "data/contact_matrices/bs1000_ngrps9_cap50_nwks5_sr33_scms.qs", 
                    "data/contact_matrices/bs1000_ngrps9_cap50_nwks5_sr41_scms.qs"),
  periods   = c('Lockdown 2', 'Lockdown 3'),
  swapouts  = 2)
{
  
  all_egs = data.table()
  i = 1
  

  cms1 = qs::qread(filenames[[1]])
  cms2 = qs::qread(filenames[[2]])
  
  
  levs <- unique(unlist(cut(seq(0,120),breaks, right=FALSE), use.names = FALSE))
  age_lab <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")
  
  # Get columns of age-groups to put into mapply
  eg = expand.grid(sort(levs),sort(levs))
  eg = data.table(eg)
  
  
  boolian_switch = matrix(eg[,Var1 == levs[swapouts] | Var2 == levs[swapouts]], nrow=81, ncol=1000)
  
  cms_aug = cms2 * data.table(boolian_switch) +  cms1 * (-(data.table(boolian_switch) -1))
  
  cms_aug
}



