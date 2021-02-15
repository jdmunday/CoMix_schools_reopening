library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot(font_size = 14) + theme(strip.background = element_blank()))


plot_timeseries_fit = function(params, title=NULL){
  
  if(is.null(title))
    {
    title = paste0('Relative_susceptibility: ', round(params$relsusc,2))
  }
  
  sigma = params$scalefactor
  i = 1.
  s = params$relsusc
  
  if(params$drop_july == TRUE){
    sr_1 = 15
  }else{
    sr_1 = 18
  }
  
  if(params$drop_peak == TRUE){
    sr_2 = 20
  }else{
    sr_2 = 22
  }
  
  if (params$sr_start < 18){
    first_part = params$sr_start:sr_1
  }else{
    first_part=NULL
  }
  
  second_part = 19:sr_2
  third_part = 23:(params$sr_end)
  fourth_part = (third_part[length(third_part)]+1):32
  
  sr_range = c(first_part, second_part, third_part, fourth_part)
  
  NWKS = length(sr_range)
  
  all_cms_all  = get_all_cms_vals('England', 'bs', 1000, 50, sr_range, nwks=2, usefirst = 100)
  
  start_dates = fread('data/start_dates.csv') 
  end_dates = fread('data/end_dates.csv')
  start_dates[, V1:=as.Date(V1, format='%d/%m/%y')]
  end_dates[, V1:=as.Date(V1, format='%d/%m/%y')]
  
  start_dates_an = sort(start_dates[survey_round %in% sr_range]$V1)
  end_dates_an = sort(end_dates[survey_round %in% (sr_range +1)]$V1)
  
  tranvec = c(i,i,i,1,1,1,1,1,1)
  suscvec = c(s,s,s,1,1,1,1,1,1)
  eigs <- compare_Rs_stan(all_cms_all, breaks, weeks_range = sr_range, suscvec = suscvec, tranvec = tranvec)
  
  
  Rsamps = eigs * sigma
  
  comix_rs = data.table(
    meanR = rowMeans(Rsamps), 
    start_dates = start_dates_an,
    end_dates = end_dates_an,
    R05 = as.vector(sapply(function(X){quantile(ecdf(Rsamps[X,]), 0.05)}, X=1:NWKS)),
    R50 = as.vector(sapply(function(X){quantile(ecdf(Rsamps[X,]), 0.5)}, X=1:NWKS)),
    R95 = as.vector(sapply(function(X){quantile(ecdf(Rsamps[X,]), 0.95)}, X=1:NWKS)), 
    epinowmean = mapply(function(X,Y){mean(england_rs[between(date, X, Y)]$mean)}, X=start_dates_an, Y=end_dates_an)
  )
  
  
  if(is.null(first_part)){
    
    first_shade_dates = c(start_dates[survey_round == second_part[1]]$V1, start_dates[survey_round == second_part[1]]$V1)
    
  }else{
    first_shade_dates =c(start_dates[survey_round == first_part[1]]$V1, end_dates[survey_round == first_part[length(first_part)]+ 1]$V1)
  }
  
  
  
  P = ggplot(comix_rs) + 
    
    annotate("rect", xmin=first_shade_dates[1], xmax=first_shade_dates[2], ymin=0.4, ymax=2,
             alpha = .2)+
    annotate("rect", xmin=start_dates[survey_round == second_part[1]]$V1, xmax=end_dates[survey_round == second_part[length(second_part)] + 1]$V1, ymin=0.4, ymax=2,
             alpha = .2)+
    annotate("rect", xmin=start_dates[survey_round == third_part[1]]$V1, xmax=end_dates[survey_round == third_part[length(third_part)]+1]$V1, ymin=0.4, ymax=2,
             alpha = .2)+
    geom_rect(aes(xmin=start_dates, xmax=end_dates, ymin=R05, ymax=R95)) + 
    geom_rect(data=comix_rs, aes(xmin=start_dates, xmax=end_dates, ymin=epinowmean-0.005, ymax=epinowmean+0.005), fill='red')  + 
    geom_line(data=england_rs[between(date, start_dates_an[1], lubridate::ymd('20210110'))], aes(x=date, y=mean)) + 
    geom_ribbon(data=england_rs[between(date, start_dates_an[1], lubridate::ymd('20210110'))], aes(x=date, ymin=lower_90, ymax=upper_90), fill='red', alpha=0.5)+
    scale_y_continuous(expand = expansion(0)) + 
    scale_x_date(breaks='month', date_labels = "%b")+
    ggtitle(title)+
    ylab('R estimate') 
  
  P
  
}

file_fits = list.files('outputs/fits/')

sigms = c()
suscs = c()
for(f in file_fits){
  
  fit = qs::qread(paste0("outputs/fits/", f))
  sigms = append(sigms, fit$par[['sigma']])
  suscs = append(suscs, fit$par[['s']])
  
}

param_dt = data.table(
  
  filename = file_fits, 
  relsusc = suscs, 
  scalefactor = sigms
  
)

param_dt[, sr_range := paste0(substring(file_fits, first = 1, last=2), ' - ', substring(file_fits, first = 4, last=5))]
param_dt[, sr_start := as.numeric(substring(file_fits, first = 1, last=2))]
param_dt[, sr_end := as.numeric(substring(file_fits, first = 4, last=5))]

param_dt[, drop_july := grepl('nojuly', filename)]
param_dt[, drop_peak := grepl('nopeak', filename)]


plot_list = list()

for(plt in 1:length(param_dt$filename)){
  
  plot_list[[plt]] = local({
    pp = plot_timeseries_fit(param_dt[plt,])
    print(pp)
    
  })
  
  
}

layout="
AEI
BFJ
CGK
DHL"

PS4 = plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_spacer() + 
  plot_list[[4]] + plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + 
  plot_list[[8]] + plot_list[[9]] + plot_list[[10]] + plot_spacer() + 
  plot_layout(design = layout)

ggsave('outputs/plots/FigureS4.pdf', PS4, width=15, height = 12)


p6 = plot_timeseries_fit(param_dt[9,], title='C')
p7 = plot_timeseries_fit(param_dt[1,], title='D')





# Plot heatmap of liklihoods ----------------

plot_likelihood_matrix =  function(ll_eg, title='NULL'){
  
  
  estimates = data.table(
    Estimate = c('1. Equal', '2. Davies et al', '3. ONS', '4. Viner et al'), 
    inf = c(1, 0.79, 1.0, 1.0), 
    susc = c(1, 0.49, 0.5, 0.64)
  )
  my_colors <- RColorBrewer::brewer.pal(4, "Pastel1")
  
  P = ggplot(ll_eg) + 
    geom_tile(aes(x = susceptibilty, y=infectiousness, fill=scaled_exp_ll)) + 
    scale_fill_gradient2(low = 'blue', mid='purple', high='orange', midpoint=0,limits=c(min(ll_eg$scaled_exp_ll), max(ll_eg$scaled_exp_ll))) + 
    geom_point(data=estimates, aes(y=inf, x=susc, shape=Estimate, color=Estimate), size=5, stroke=2)+
    scale_color_manual(values = my_colors, name='Study:', guide=FALSE)+
    scale_shape_discrete(name = "Study:")+
    guides(fill=guide_colorbar(title='exp(-ll/1000)', label.position='top', barwidth = 15, title.position = 'bottom'), shape=guide_legend(ncol = 4, label.position = 'top'))+
    xlab('Relative Susceptibility') + 
    ylab('Relative Infectiousness') +
    ggtitle(title)+ 
    theme(legend.position = "top", legend.direction = "horizontal", legend.box = 'vertical', legend.box.just = 'left', legend.justification = 'left')
  
  P
  
}

ll_eg = qs::qread('outputs/fitmats/22_19_28_llmat.qs')
p3 = plot_likelihood_matrix(ll_eg, title='A')

ll_eg = qs::qread('outputs/fitmats/26_12_28_llmat.qs')
p4 = plot_likelihood_matrix(ll_eg, title='B')



layout <- "
DDDDDDEEEEEE
DDDDDDEEEEEE
DDDDDDEEEEEE
FFFFFFFFFFFF
FFFFFFFFFFFF
GGGGGGGGGGGG
GGGGGGGGGGGG
"

Fig1 = p3 + p4 + p6 + p7 +
  plot_layout(design = layout)

ggsave('outputs/plots//Figure1.pdf', Fig1, width = 10, height = 14)





layout <- "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
"

Fig2 = p1 + p2 + 
  plot_layout(design = layout)

ggsave('outputs/plots/Figure2.pdf', Fig2, width = 15, height = 10)
