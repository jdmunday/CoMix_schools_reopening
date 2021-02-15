library(wesanderson)
all_egs = get_all_egs()


all_egs[,cms_lab := round(cms, 1)]

switchout = 3

cms_aug = augment_cms(swapout=switchout)

levs <- unique(unlist(cut(seq(0,120),breaks, right=FALSE), use.names = FALSE))
age_lab <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")

# Get columns of age-groups to put into mapply
eg = expand.grid(sort(levs),sort(levs))
eg = data.table(eg)
eg[, Var1 := factor(Var1, levels = levs, labels = age_lab)]
eg[, Var2 := factor(Var2, levels = levs, labels = age_lab)]

eg[,cms := rowMeans(cms_aug)]

eg[,sr:='Open Primary']
#eg[, daterange := paste0(dts$start_date, ' - ', dts$end_date)]

eg[,cms_lab := round(cms, 1)]

all_egs = rbind(all_egs, eg)

all_cms_NS = list()

all_cms_NS[[1]] = as.matrix(cms_aug)




switchout = 2

cms_aug = augment_cms(swapout=switchout)

levs <- unique(unlist(cut(seq(0,120),breaks, right=FALSE), use.names = FALSE))
age_lab <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")

# Get columns of age-groups to put into mapply
eg = expand.grid(sort(levs),sort(levs))
eg = data.table(eg)
eg[, Var1 := factor(Var1, levels = levs, labels = age_lab)]
eg[, Var2 := factor(Var2, levels = levs, labels = age_lab)]

eg[,cms := rowMeans(cms_aug)]

eg[,sr:='Open Secondary']
#eg[, daterange := paste0(dts$start_date, ' - ', dts$end_date)]

eg[,cms_lab := round(cms, 1)]

all_egs = rbind(all_egs, eg)

all_egs[sr == 'Lockdown 2', sr:='Open Both']
theme_set(cowplot::theme_cowplot(font_size = 24) + theme(strip.background = element_blank()))

ggplot(all_egs[sr %in% list('Open Primary', 'Open Secondary','Open Both')], aes(Var1, Var2, fill= cms, label=cms_lab)) + 
  geom_tile()+
  geom_text(color='white', size=5)+
  facet_wrap( ~sr, ncol=1) +
  scale_fill_viridis(discrete=FALSE, name='Mean \ncontacts', begin=0, end=1., limits = c(0,4.0))+ 
  ylab('Contact age group') +
  xlab('Participant age group') +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank())


library(ggridges)

all_cms_NP = list()

all_cms_NP[[1]] = as.matrix(cms_aug)
eigs_LK2NS = get_r_estimates('England', 'bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NS)


eigs_LK2NP = get_r_estimates('England', 'bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NP)


eigs_lK3AS =  get_r_estimates('England', 'bs', 1000, 50, 41:41, nwks=3)
eigs_lK2AS = get_r_estimates('England', 'bs', 1000, 50, 33:33, nwks=5)


eig_stats = data.table()
periods = c('Only Primary', 'Only Secondary', 'All Schools', 'No Schools')
suscrule = rep('Equal', 4)
i = 1
for(rats in list(eigs_LK2NP, eigs_LK2NP, eigs_lK2AS, eigs_lK3AS)){
  eig_stats = rbind(eig_stats, data.table(
    
    compare = periods[i], 
    suscrule = suscrule[i],
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


eig_rats = data.table()
periods = c('Only Primary', 'Only Secondary', 'All Schools')
suscrule = rep('Equal', 3)
i = 1
for(rats in list(l2NSRAT, l2NPRAT, l2ASRAT)){
  eig_rats = rbind(eig_rats, data.table(
    
    compare = periods[i], 
    suscrule = suscrule[i],
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



suscvec_ = c(0.4,
            0.4,
            0.4,
            0.79,
            0.86,
            0.8,
            0.82,
            0.88,
            0.74)

tranvec_ = c(0.645,
            0.645,
            0.605,
            0.635,
            0.665,
            0.7,
            0.745,
            0.815,
            0.845)


eigs_LK2NS = get_r_estimates('England', 'bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NS, suscvec = suscvec_, tranvec = tranvec_)

eigs_LK2NP = get_r_estimates('England', 'bs', 1000, 50, 1:1, nwks=8, cms_all=all_cms_NP, suscvec = suscvec_, tranvec = tranvec_)


eigs_lK3AS =  get_r_estimates('England', 'bs', 1000, 50, 41:41, nwks=3, suscvec = suscvec_, tranvec = tranvec_)
eigs_lK2AS = get_r_estimates('England', 'bs', 1000, 50, 33:33, nwks=5, suscvec = suscvec_, tranvec = tranvec_)

eig_stats_ST = data.table()
periods = c('Only Primary', 'Only Secondary', 'All Schools', 'No Schools')
suscrule = rep('Davies et al', 4)
i = 1
for(rats in list(eigs_LK2NP, eigs_LK2NP, eigs_lK2AS, eigs_lK3AS)){
  eig_stats_ST = rbind(eig_stats_ST, data.table(
    
    compare = periods[i], 
    suscrule = suscrule[i],
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

l2NPRAT_ST = eigs_LK2NP / eigs_lK3AS
l2NSRAT_ST = eigs_LK2NS / eigs_lK3AS
l2ASRAT_ST = eigs_lK2AS / eigs_lK3AS



eig_rats_ST = data.table()
periods = c('Only Primary', 'Only Secondary', 'All Schools')
suscrule = rep('Davies et al', 3)
i = 1
for(rats in list(l2NSRAT_ST, l2NPRAT_ST, l2ASRAT_ST)){
  eig_rats_ST = rbind(eig_rats_ST, data.table(
    
    compare = periods[i], 
    suscrule = suscrule[i],
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


fmt <- function(x) {
  format(round(x, 1), nsmall = 1)
}

write.csv(eig_stats, '../comix/lockdown_eigs_plain.csv')
write.csv(eig_rats, '../comix/lockdown_rats_plain.csv')
write.csv(eig_stats_ST, '../comix/lockdown_eigs_st.csv')
write.csv(eig_rats_ST, '../comix/lockdown_rats_st.csv')

suscrule = c('Equal','Equal','Equal', 'Davies', 'Davies', 'Davies')


fmt <- function(x) {
  format(round(x, 1), nsmall = 1)
}
rat_table = rbind(eig_rats, eig_rats_ST)

rat_table[, ratio := paste0(fmt(rat50), ' (', fmt(rat05), ' - ', fmt(rat95), ')' )]

slimtab = rat_table[,c('suscrule', 'compare', 'ratio')]

colnames(slimtab) = c('Susceptibility/ Infectiousness', 'Attendence', 'Ratio')






theme_set(cowplot::theme_cowplot(font_size = 20) + theme(strip.background = element_blank()), pal)

eig_frame = data.table()
periods = c('Only Primary - Equal', 'Only Secondary - Equal', 'All Schools - Equal', 'Only Primary - D', 'Only Secondary - D', 'All Schools - D')
suscrule = c('Equal','Equal','Equal', 'Davies', 'Davies', 'Davies')
i = 1
for(rats in list(l2NSRAT, l2NPRAT, l2ASRAT, l2NSRAT_ST, l2NPRAT_ST, l2ASRAT_ST)){
  eig_frame = rbind(eig_frame, data.table(
    
    Period = periods[i], 
    Eigenvalue = as.vector(rats),
    suscrule = suscrule[i]
    
  ))
  i = i + 1
}

pal <- wes.palette(name = "Zissou", type = "continuous")

ggplot(data = eig_frame[Period %in% periods], aes(x=Eigenvalue, y=Period)) + 
  geom_density_ridges(aes(fill=suscrule), color='white' ) + 
  scale_fill_brewer(palette="Set1")





epinow = fread('../comix/data/epinowouts.csv')
