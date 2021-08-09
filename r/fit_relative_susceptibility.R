library(data.table)
source('r/get_contact_matrices.R')
source('r/compare_Rs.R')

epinow = data.table::fread('data/rt_estimates/rt_cases.csv')
england_rs = epinow[region=='England']
england_rs[, date := as.Date(date)]

breaks=c(0,5,12,18,30,40,50,60,70,Inf)

opt_func_multi = function(x, par, i){
  
  s = par[['s']]
  sigma = par[['sigma']]
  tranvec = c(i,i,i,1,1,1,1,1,1)
  suscvec = c(s,s,s,1,1,1,1,1,1)
  eigs <- compare_Rs_stan(all_cms_all, breaks, weeks_range = sr_range, suscvec = suscvec, tranvec = tranvec)
  
  
  Ms = mapply(function(X,Y){mean(england_rs[between(date, X, Y)]$mean)}, X=start_dates_an, Y=end_dates_an)
  Ss = mapply(function(X,Y){mean(england_rs[between(date, X, Y)]$sd)}, X=start_dates_an, Y=end_dates_an)
  
  gamma_params = mapply(function(X,Y){ConnMatTools::gammaParamsConvert(mean=X, sd=Y)}, X=Ms, Y=Ss)
  
  
  Rsamps = eigs * sigma
  
  
  ll = mapply(function(X,Y){dgamma(x=Y, shape=gamma_params[,X]$shape, scale=gamma_params[,X]$scale, log = TRUE)}, X=1:NWKS, Y=Rsamps)
  
  
  -(sum(ll))
  
}

bbm_func_multi = function(s, sigma, i){
  
  #s = par[['s']]
  #sigma = par[['sigma']]
  tranvec = c(i,i,i,1,1,1,1,1,1)
  suscvec = c(s,s,s,1,1,1,1,1,1)
  eigs <- compare_Rs_stan(all_cms_all, breaks, weeks_range = sr_range, suscvec = suscvec, tranvec = tranvec)
  
  
  Ms = mapply(function(X,Y){mean(england_rs[between(date, X, Y)]$mean)}, X=start_dates_an, Y=end_dates_an)
  Ss = mapply(function(X,Y){mean(england_rs[between(date, X, Y)]$sd)}, X=start_dates_an, Y=end_dates_an)
  
  gamma_params = mapply(function(X,Y){ConnMatTools::gammaParamsConvert(mean=X, sd=Y)}, X=Ms, Y=Ss)
  
  
  Rsamps = eigs * sigma
  
  
  ll = mapply(function(X,Y){dgamma(x=Y, shape=gamma_params[,X]$shape, scale=gamma_params[,X]$scale, log = TRUE)}, X=1:NWKS, Y=Rsamps)
  
  
  -(sum(ll))
  
}




sr_range = c(12:15, 19:20, 23:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)

start_dates = fread('data/start_dates.csv') 
end_dates = fread('data/end_dates.csv')
start_dates[, V1:=as.Date(V1, format='%d/%m/%y')]
end_dates[, V1:=as.Date(V1, format='%d/%m/%y')]


start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)

outs10_multi = optim(c(s = 0.3, sigma=0.5), lower=c(1e-5, 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)
bbouts_10 = bbmle::mle2(bbm_func_multi, start =list(s = 0.3, sigma=0.5), fixed = list(i=1.0))

summary(bbouts_10)
bbprof_10 = bbmle::profile(bbouts_10)

qs::qsave(outs10_multi,'outputs/fits/12_28_nojuly_nopeak.qs')






sr_range = c(12:15, 19:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)

outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)
bbouts_10 = bbmle::mle2(bbm_func_multi, start =list(s = 0.3, sigma=0.5), fixed = list(i=1.0))

bbmle::summary(bbouts_10)
bbprof_10 = bbmle::profile(bbouts_10)
bbprof_10
qs::qsave(outs10_multi,'outputs/fits/12_28_nojuly.qs')

sr_range = c(12:20, 23:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)

outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/12_28_nopeak.qs')





sr_range = c(12:33)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/12_33_alldata.qs')



sr_range = c(12:15, 19:20, 23:33)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/12_33_nojuly_nopeak.qs')



sr_range = c(12:15, 19:33)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/12_33_nojuly.qs')


sr_range = c(12:20, 23:33)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/12_33_nopeak.qs')


sr_range = c(19:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/19_28_alldata.qs')

sr_range = c(19:20,23:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)



start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)
bbouts_10 = bbmle::mle2(bbm_func_multi, start =list(s = 0.3, sigma=0.5), fixed = list(i=1.0))

bbprof_10 = bbmle::confint(bbouts_10)

qs::qsave(outs10_multi,'outputs/fits/19_28_nopeak.qs')

sr_range = c(19:20,23:33)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)


start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
outs10_multi = optim(c(s = 0.3, sigma=0.5), lower = c(s = 1e-5, sigma = 1e-5), fn = opt_func_multi, method = "L-BFGS-B", x=c(0,0), i=1.0)

qs::qsave(outs10_multi,'outputs/fits/19_33_nopeak.qs')





sr_range = c(19:20, 23:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)


start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)

all_ll = data.table()
for (s in seq(0.3,1.2,0.05)){
  row_ll = c()
  for (i in seq(0.3,1.2,0.05)){
    row_ll = append(row_ll,opt_func_multi(c(0,0), par=c(s=s, sigma=0.22), i=i))
  }
  all_ll = rbind(all_ll, row_ll)   
}


ll_eg= expand.grid(seq(0.3,1.2,0.05), seq(0.3,1.2,0.05))
colnames(ll_eg) = c('susceptibilty', 'infectiousness')
ll_eg = data.table(ll_eg)
ll_eg[, scaled_exp_ll := exp(-all_ll$x/1000)]
qs::qsave(ll_eg,'outputs/fitmats/22_19_28_llmat.qs')


sr_range = c(12:15, 19:20, 23:28)
NWKS = length(sr_range)

all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)


start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)

all_ll = data.table()
for (s in seq(0.3,1.2,0.05)){
  row_ll = c()
  for (i in seq(0.3,1.2,0.05)){
    row_ll = append(row_ll,opt_func_multi(c(0,0), par=c(s=s, sigma=0.26), i=i))
  }
  all_ll = rbind(all_ll, row_ll)   
}


ll_eg= expand.grid(seq(0.3,1.2,0.05), seq(0.3,1.2,0.05))
colnames(ll_eg) = c('susceptibilty', 'infectiousness')
ll_eg = data.table(ll_eg)
ll_eg[, scaled_exp_ll := exp(-all_ll$x/1000)]
qs::qsave(ll_eg,'outputs/fitmats/26_12_28_llmat.qs')






