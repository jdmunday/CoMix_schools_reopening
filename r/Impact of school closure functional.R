library(viridis)
library(ggplot2)
library(wesanderson)
library(data.table)
source('r/get_contact_matrices.R')
source('r/get_eigenvalues.R')
source('r/compare_Rs.R')



breaks = c(0,5,12,18,30,40,50,60,70,Inf)


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


all_cms_NP = list()

all_cms_NP[[1]] = as.matrix(cms_aug)


suscvec_davies = c(0.4,
             0.4,
             0.4,
             0.79,
             0.86,
             0.8,
             0.82,
             0.88,
             0.74)

tranvec_davies = c(0.645,
             0.645,
             0.605,
             0.635,
             0.665,
             0.7,
             0.745,
             0.815,
             0.845)


suscvec_ons2= c(0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0)

tranvec_ons2 = c(1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0)

suscvec_meta = c(0.64 ,
                 0.64 ,
                 0.64 ,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0)

tranvec_meta = c(1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0)

suscvec_fit = c(0.31 ,
                 0.31 ,
                 0.31 ,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0)

tranvec_fit = c(1.0,
                1.0,
                1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0)




eig_stats = data.table()
eig_rats = data.table()
eig_frame = data.table()

eig_frames = list(eig_stats, eig_rats, eig_frame)

eig_frames = calculate_eigs_and_rats(eig_stats_ = eig_frames[[1]], 
                                     eig_rats_ = eig_frames[[2]],
                                     eig_frame_ = eig_frames[[3]], 
                                     all_cms_NP = all_cms_NP,
                                     all_cms_NS = all_cms_NS)

eig_frames = calculate_eigs_and_rats(eig_stats_ = eig_frames[[1]], 
                                     eig_rats_ = eig_frames[[2]],  
                                     eig_frame_ = eig_frames[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '2. Davies et al',
                                     suscvec_ = suscvec_davies,
                                     tranvec_ = tranvec_davies)




eig_frames = calculate_eigs_and_rats(eig_stats_ = eig_frames[[1]], 
                                     eig_rats_ = eig_frames[[2]],  
                                     eig_frame_ = eig_frames[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '3. ONS', 
                                     suscvec_ = suscvec_ons2, 
                                     tranvec_ = tranvec_ons2)

eig_frames = calculate_eigs_and_rats(eig_stats_ = eig_frames[[1]], 
                                     eig_rats_ = eig_frames[[2]],  
                                     eig_frame_ = eig_frames[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '4. Viner et al', 
                                     suscvec_ = suscvec_meta, 
                                     tranvec_ = tranvec_meta)

eig_frames = calculate_eigs_and_rats(eig_stats_ = eig_frames[[1]], 
                                     eig_rats_ = eig_frames[[2]],  
                                     eig_frame_ = eig_frames[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '5. CoMix fit', 
                                     suscvec_ = suscvec_fit, 
                                     tranvec_ = tranvec_fit)

theme_set(cowplot::theme_cowplot(font_size = 14) + theme(strip.background = element_blank()))
library(ggridges)
p1 = ggplot(data = eig_frames[[3]], aes(x=Ratio, y=suscrule)) + 
  facet_wrap(~ Period, ncol=1)+
  geom_vline(xintercept = 1., size=1.5, linetype='dashed')+
  geom_density_ridges(aes(fill=suscrule), color='white') + 
  ylab("") + 
  xlab("Relative increase in R")+
  guides(fill=guide_legend(label.position = 'top', ))+
  scale_y_discrete( labels=c('1.','2.','3.','4.', '5.'))+
  scale_fill_brewer(palette="Pastel1", name = "Study:", type = 'cols')+
  ggtitle('A') + 
  theme(
    legend.position = "top"
    
    ) 
fmt = function(x) {
  format(round(x, 1), nsmall = 1)
}

rat_table = eig_frames[[2]]

rat_table[, ratio := paste0(fmt(rat50), ' (', fmt(rat05), ' - ', fmt(rat95), ')' )]

slimtab = rat_table[,c('suscrule', 'compare', 'ratio')]

colnames(slimtab) = c('SuscInf', 'Attendence', 'Ratio')

slimtab

write.csv(slimtab, 'outputs/tables//eigenvalue_ratios.csv')

R_table = data.table()
for (R in c(0.7, 0.8, 0.9, 1.0)){
  print(R)
  new_rs = R * rat_table[,c('rat05',    'rat50',    'rat95')]
  colnames(new_rs) = c('newR05',    'R0',    'newR95')
  new_rs[, Rold:= R]
  new_rs[, Roldtag:= paste0('Baseline R = ',R)]
  r_table = cbind(slimtab,new_rs)
  R_table = rbind(R_table, r_table)
}

R_table

R_table_pres = data.table(
  Scenario = R_table$SuscInf, 
  Attendance = R_table$Attendence, 
  Baseline_R = R_table$Rold, 
  New_R = paste0(fmt(R_table$R0), ' (', fmt(R_table$newR05), ' - ', fmt(R_table$newR95), ')')
)

R_table_pres = dcast(R_table_pres, Scenario + Attendance ~ Baseline_R, value.var = 'New_R')

write.csv(R_table_pres, 'outputs/tables/R_estimates_new.csv')

p2 = ggplot(R_table) + 
  geom_pointrange(mapping=aes(y=Attendence, x=R0, xmin=newR05, xmax=newR95, color=SuscInf, group=SuscInf), position = position_dodge(width=0.5),size=0.5) + 
  #geom_point(mapping=aes(y=Attendence, x=Rold, color=SuscInf, group=SuscInf),position = position_dodge(width=1),  size=2.5) + 
  geom_vline(aes(xintercept = Rold), color='MidnightBlue', alpha=0.5)+
  geom_vline(xintercept = 1.0, color='darkgrey', linetype = "dashed", alpha=0.8)+
  facet_wrap(~ Roldtag, ncol=1) + 
  guides(color=guide_legend(title="Susceptibility and infectiousness profile", nrow=2)) + 
  scale_color_brewer(palette="Pastel1") + 
  scale_x_continuous(breaks = seq(0.8, 4.0, 0.2)) + 
  ylab('')+
  xlab('R')+
  ggtitle('B')+ 
  scale_y_discrete(limits = unique((R_table$Attendence)))+
  theme(legend.position = "none")
  
  
library(patchwork)

layout <- "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
"

Fig2 = p1 + p2 + 
  plot_layout(design = layout)

ggsave('outputs/plots/Figure2.png', Fig2, width = 15, height = 10)


theme_set(cowplot::theme_cowplot(font_size = 20) + theme(strip.background = element_blank()))

BP = ggplot(R_table[Rold %in% c(0.8)]) + 
  geom_pointrange(mapping=aes(y=Attendence, x=R0, xmin=newR05, xmax=newR95, color=SuscInf, group=SuscInf), position = position_dodge(width=0.5),size=1) + 
  #geom_point(mapping=aes(y=Attendence, x=Rold, color=SuscInf, group=SuscInf),position = position_dodge(width=1),  size=2.5) + 
  geom_vline(aes(xintercept = Rold), color='MidnightBlue', alpha=0.5)+
  geom_vline(xintercept = 1.0, color='darkgrey', linetype = "dashed", alpha=0.8)+ 
  guides(color=guide_legend(title="Susceptibility/ \ninfectiousness \nestimate", nrow=2)) + 
  scale_colour_manual(values = wes_palette("Darjeeling1")) + 
  scale_x_continuous(breaks = seq(0.8, 4.0, 0.2)) + 
  ylab('Reopening strategy')+
  xlab('R')+
  ggtitle('Potential impact on R of opening schools \n(Assuming current R is 0.8 with schools closed)') + 
  scale_y_discrete(limits = unique((R_table$Attendence)), labels=c('Primary', 'Secondary', 'All Schools'))+
  theme(legend.position = "bottom")
ggsave('outputs/plots/blogplot_08.png', BP, width=10, height=8)


PP = ggplot(R_table) + 
  geom_pointrange(mapping=aes(y=Attendence, x=R0, xmin=newR05, xmax=newR95, color=SuscInf, group=SuscInf), position = position_dodge(width=0.5),size=1) + 
  #geom_point(mapping=aes(y=Attendence, x=Rold, color=SuscInf, group=SuscInf),position = position_dodge(width=1),  size=2.5) + 
  geom_vline(aes(xintercept = Rold), color='MidnightBlue', alpha=0.5)+
  geom_vline(xintercept = 1.0, color='darkgrey', linetype = "dashed", alpha=0.8)+
  facet_wrap(~ Roldtag, ncol=1) +  
  guides(color=guide_legend(title="Susceptibility/ \ninfectiousness \nestimate", nrow=2)) + 
  scale_colour_manual(values = wes_palette("Darjeeling1")) + 
  scale_x_continuous(breaks = seq(0.8, 4.0, 0.2)) + 
  ylab('Reopening strategy')+
  xlab('R')+
  ggtitle('Potential impact on R of opening schools') + 
  scale_y_discrete(limits = unique((R_table$Attendence)), labels=c('Primary', 'Secondary', 'All Schools'))+
  theme(legend.position = "bottom")
ggsave('outputs/plots/presplot.png', PP, width=10, height=10)


PP1 = ggplot(data = eig_frames[[3]], aes(x=Ratio, y=suscrule)) + 
  facet_wrap(~ Period, ncol=1)+
  geom_vline(xintercept = 1., size=1.5, linetype='dashed')+
  geom_density_ridges(aes(fill=suscrule), color='white') + 
  ylab("") + 
  xlab("Relative increase in R")+
  guides(fill=guide_legend(title="Susceptibility/ \ninfectiousness \nestimate", nrow=2))+
  scale_y_discrete( labels=c('1.','2.','3.','4.', '5.'))+
  scale_fill_manual(values = wes_palette("Darjeeling1")) + 
  ggtitle('Relative increase in R upon reopening schools') + 
  theme(
    legend.position = "bottom")
ggsave('outputs/plots/presplot1.png', PP1, width=10, height=10)






#------------

vax1 = c(1.0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.5, 0.2)
vax2 = c(1.0, 1.0, 1.0, 0.9, 0.9, 0.5, 0.2, 0.2, 0.1)
vax3 = c(1.0, 1.0, 1.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1)

eig_stats = data.table()
eig_rats = data.table()
eig_frame = data.table()

eig_frames_vax = list(eig_stats, eig_rats, eig_frame)


eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                     eig_rats_ = eig_frames_vax[[2]],  
                                     eig_frame_ = eig_frames_vax[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '2. Davies et al - vax 1',
                                     suscvec_ = suscvec_davies * vax1,
                                     tranvec_ = tranvec_davies)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '2. Davies et al - vax 2',
                                         suscvec_ = suscvec_davies  * vax2,
                                         tranvec_ = tranvec_davies)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '2. Davies et al - vax 3',
                                         suscvec_ = suscvec_davies  * vax3,
                                         tranvec_ = tranvec_davies)


eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                     eig_rats_ = eig_frames_vax[[2]],  
                                     eig_frame_ = eig_frames_vax[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '3. ONS - vax 1', 
                                     suscvec_ = suscvec_ons2 * vax1, 
                                     tranvec_ = tranvec_ons2)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '3. ONS - vax 2', 
                                         suscvec_ = suscvec_ons2 * vax2, 
                                         tranvec_ = tranvec_ons2)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '3. ONS - vax 3', 
                                         suscvec_ = suscvec_ons2 * vax3, 
                                         tranvec_ = tranvec_ons2)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                     eig_rats_ = eig_frames_vax[[2]],  
                                     eig_frame_ = eig_frames_vax[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '4. Viner et al - vax 1', 
                                     suscvec_ = suscvec_meta * vax1, 
                                     tranvec_ = tranvec_meta)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '4. Viner et al - vax 2', 
                                         suscvec_ = suscvec_meta * vax2, 
                                         tranvec_ = tranvec_meta)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '4. Viner et al - vax 3', 
                                         suscvec_ = suscvec_meta * vax3, 
                                         tranvec_ = tranvec_meta)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                     eig_rats_ = eig_frames_vax[[2]],  
                                     eig_frame_ = eig_frames_vax[[3]],  
                                     all_cms_NP = all_cms_NP, 
                                     all_cms_NS = all_cms_NS, 
                                     suscname = '5. CoMix fit - vax 1', 
                                     suscvec_ = suscvec_fit * vax1, 
                                     tranvec_ = tranvec_fit)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '5. CoMix fit - vax 2', 
                                         suscvec_ = suscvec_fit * vax2, 
                                         tranvec_ = tranvec_fit)

eig_frames_vax = calculate_eigs_and_rats(eig_stats_ = eig_frames_vax[[1]], 
                                         eig_rats_ = eig_frames_vax[[2]],  
                                         eig_frame_ = eig_frames_vax[[3]],  
                                         all_cms_NP = all_cms_NP, 
                                         all_cms_NS = all_cms_NS, 
                                         suscname = '5. CoMix fit - vax 3', 
                                         suscvec_ = suscvec_fit * vax3, 
                                         tranvec_ = tranvec_fit)

pvax = ggplot(data = eig_frames_vax[[3]], aes(x=Ratio, y=suscrule)) + 
  facet_wrap(~ Period, ncol=1)+
  geom_vline(xintercept = 1., size=1.5, linetype='dashed')+
  geom_density_ridges(aes(fill=suscrule), color='white') + 
  ylab("") + 
  xlab("Relative increase in R")+
  guides(fill=guide_legend(label.position = 'top'))+
  scale_color_discrete(name='Scenario')+
  scale_y_discrete( labels=c())+
  ggtitle('A')
  
  
vax = data.table(
    age = eg$Var1[1:9],
    vax1 = 1 - vax1+0.01, 
    vax2 = 1 - vax2+0.01, 
    vax3 = 1 - vax3+0.01
  )

vax_long = 
  melt(vax, id.vars = ('age'), 
     measure.vars = c('vax1', 
                      'vax2', 
                      'vax3'), 
     value.name = 'vax', 
     variable.name = 'scenario')

vax_bar = ggplot(vax_long, aes(x=age, y=vax)) + geom_bar(stat="identity") + facet_wrap('scenario', ncol =1)
                 
                 
eigs_no_vax_davies = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_davies, tranvec = tranvec_davies)
eigs_vax_1_davies = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_davies * vax1, tranvec = tranvec_davies)
eigs_vax_2_davies = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_davies * vax2, tranvec = tranvec_davies)
eigs_vax_3_davies = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_davies * vax3, tranvec = tranvec_davies)

relreduct_vax1_davies  = eigs_vax_1_davies  / eigs_no_vax_davies 
relreduct_vax2_davies  = eigs_vax_2_davies  / eigs_no_vax_davies 
relreduct_vax3_davies  = eigs_vax_3_davies  / eigs_no_vax_davies 

eigs_no_vax_ONS = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5, suscvec = suscvec_ons2,        tranvec = tranvec_ons2)
eigs_vax_1_ONS = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_ons2 * vax1, tranvec = tranvec_ons2)
eigs_vax_2_ONS = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_ons2 * vax2, tranvec = tranvec_ons2)
eigs_vax_3_ONS = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_ons2 * vax3, tranvec = tranvec_ons2)

relreduct_vax1_ONS = eigs_vax_1_ONS / eigs_no_vax_ONS
relreduct_vax2_ONS = eigs_vax_2_ONS / eigs_no_vax_ONS
relreduct_vax3_ONS = eigs_vax_3_ONS / eigs_no_vax_ONS

eigs_no_vax_viner = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5, suscvec = suscvec_meta,        tranvec = tranvec_meta)
eigs_vax_1_viner = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_meta * vax1, tranvec = tranvec_meta)
eigs_vax_2_viner = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_meta * vax2, tranvec = tranvec_meta)
eigs_vax_3_viner = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_meta * vax3, tranvec = tranvec_meta)

relreduct_vax1_viner = eigs_vax_1_viner / eigs_no_vax_viner
relreduct_vax2_viner = eigs_vax_2_viner / eigs_no_vax_viner
relreduct_vax3_viner = eigs_vax_3_viner / eigs_no_vax_viner

eigs_no_vax_fit = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5, suscvec = suscvec_fit,        tranvec = tranvec_fit)
eigs_vax_1_fit = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_fit * vax1, tranvec = tranvec_fit)
eigs_vax_2_fit = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_fit * vax2, tranvec = tranvec_fit)
eigs_vax_3_fit = get_eigenvalues('bs', 1000, 50, 41:41, nwks=5,  suscvec = suscvec_fit * vax3, tranvec = tranvec_fit)

relreduct_vax1_fit = eigs_vax_1_fit / eigs_no_vax_fit
relreduct_vax2_fit = eigs_vax_2_fit / eigs_no_vax_fit
relreduct_vax3_fit = eigs_vax_3_fit / eigs_no_vax_fit
i=1
suscs = c(rep(c('2. Davies et al'), 3),rep(c('3. ONS' ), 3), rep(c('4. Viner et al'), 3), rep(c('5. CoMix fit'), 3))
vaxs = rep(c('Vacc 1', 'Vacc 2', 'Vacc 3'), 4)
prop_rats = data.table()
for(rats in list(relreduct_vax1_davies, 
                 relreduct_vax2_davies, 
                 relreduct_vax3_davies,
                 relreduct_vax1_ONS,
                 relreduct_vax2_ONS,
                 relreduct_vax3_ONS,
                 relreduct_vax1_viner,
                 relreduct_vax2_viner,
                 relreduct_vax3_viner,
                 relreduct_vax1_fit,
                 relreduct_vax2_fit,
                 relreduct_vax3_fit
                 )){
  prop_rats = rbind(prop_rats, data.table(
    vaccination = vaxs[i], 
    suscrule = suscs[i],
    compare_num = i,
    rat05 = quantile(ecdf(rats), 0.05),
    rat50 = quantile(ecdf(rats), 0.5),
    rat95 = quantile(ecdf(rats), 0.95),
    rat25 = quantile(ecdf(rats), 0.25),
    rat75 = quantile(ecdf(rats), 0.75),
    rat10 = quantile(ecdf(rats), 0.10),
    rat90 = quantile(ecdf(rats), 0.9)
  ))
  i = i+1
}
eig_frames_vax_2 = data.table(eig_frames_vax[[2]])

str_front_trunc= function(x,n){substr(x, start = 1+n, stop=nchar(x))}
str_end_trunc= function(x,n){substr(x, start = 1, stop=nchar(x)-n)}

split_n_trunc_f = function(x){str_front_trunc(strsplit(x, '-')[[1]][2], 1)}
split_n_trunc_e = function(x){str_end_trunc(strsplit(x, '-')[[1]][1], 1)}
suslist = eig_frames_vax_2$suscrule

eig_frames_vax_2[, vaccination := sapply(suslist, split_n_trunc_f)]
eig_frames_vax_2[, suscrule := sapply(suslist, split_n_trunc_e)]


eig_frames_vax_2[, ratio := paste0(fmt(rat50), ' (', fmt(rat05), ' - ', fmt(rat95), ')' )]

slimtab_vax = eig_frames_vax_2[,c('suscrule', 'vaccination', 'compare', 'ratio')]

colnames(slimtab_vax) = c('suscrule', 'vaccination', 'Attendence', 'Ratio')

prop_rats[, vaxred := paste0(fmt(rat50), ' (', fmt(rat05), ' - ', fmt(rat95), ')' )]
slim_prop_rats = prop_rats[,c('suscrule', 'vaccination', 'vaxred')]
colnames(slim_prop_rats) = c('suscrule', 'vaccination',  'reduction')

merge(slimtab_vax[Attendence == 'Both'], slim_prop_rats, by = c('suscrule', 'vaccination'))


prop_frame = data.table()
i = 1
for(rats in list(relreduct_vax1_davies, 
                 relreduct_vax2_davies, 
                 relreduct_vax3_davies,
                 relreduct_vax1_ONS,
                 relreduct_vax2_ONS,
                 relreduct_vax3_ONS,
                 relreduct_vax1_viner,
                 relreduct_vax2_viner,
                 relreduct_vax3_viner,
                 relreduct_vax1_fit,
                 relreduct_vax2_fit,
                 relreduct_vax3_fit
)){
  prop_frame = rbind(prop_frame, data.table(
    
    vaccination = rep(vaxs[i], length(rats)), 
    Ratio = as.vector(rats),
    suscrule = rep(suscs[i], length(rats))
    
  ))
  i = i + 1
}

prop_frame[, School_Opening := eig_frames_vax_3[Period=='Both']$Ratio]
prop_frame[, Overall_impact := Ratio * School_Opening]

oi_summary = data.table()

for (susc in suscs[seq(1,12,3)]){
  for (vax in vaxs[1:3]){
    ecdf_rats = ecdf(prop_frame[vaccination==vax & suscrule==susc]$Overall_impact)
    oi_summary = rbind(oi_summary, 
                         data.table(
                           vaccination = vax, 
                           suscrule = susc,
                           rat05 = quantile(ecdf_rats, 0.05),
                           rat50 = quantile(ecdf_rats, 0.5),
                           rat95 = quantile(ecdf_rats, 0.95),
                           rat25 = quantile(ecdf_rats, 0.25),
                           rat75 = quantile(ecdf_rats, 0.75),
                           rat10 = quantile(ecdf_rats, 0.10),
                           rat90 = quantile(ecdf_rats, 0.9)
                           
                           
                         ))
    
    
  }}


oi_summary[, overallimpact := paste0(fmt(rat50), ' (', fmt(rat05), ' - ', fmt(rat95), ')' )]
oi_summary_slim = oi_summary[,c('suscrule', 'vaccination', 'overallimpact')]

impact_tab = merge(slimtab_vax[Attendence == 'Both'], slim_prop_rats, by = c('suscrule', 'vaccination'))

impact_tab = merge(impact_tab, oi_summary_slim,  by = c('suscrule', 'vaccination'))
names(impact_tab) = c('Study', 'Vaccination scenario', 'Attendance', 'Rel. Inc. - Schools', 'Rel. Red. - Vaccination', 'Overall Impact')

write.csv(impact_tab, 'vacc_schools_combined_impact.csv')

pvax1 = ggplot(data = prop_frame, aes(x=School_Opening, y=suscrule)) + 
  facet_wrap(~ vaccination, ncol=1)+
  geom_vline(xintercept = 1., size=1.5, linetype='dashed')+
  geom_density_ridges(aes(fill=suscrule), color='white') + 
  ylab("") + 
  xlab("Relative increase in R")+
  ggtitle('A')+
  theme(
    legend.position = "None"
  )

pvax2 = ggplot(data = prop_frame, aes(x=Overall_impact, y=suscrule)) + 
  facet_wrap(~ vaccination, ncol=1)+
  geom_vline(xintercept = 1., size=1.5, linetype='dashed')+
  geom_density_ridges(aes(fill=suscrule), color='white') + 
  ylab("") + 
  xlab("Relative increase in R")+
  guides(fill=guide_legend(label.position = 'top'))+
  scale_color_discrete(name='Scenario')+
  scale_y_discrete( labels=c())+
  ggtitle('B')+
  theme(
    legend.title = element_blank()
  )


pvax1 + pvax2

rat_table_mod = copy(rat_table)
rat_table_mod[, overallimpact := ratio]
rat_table_mod = rat_table_mod[compare == 'Both']
rat_table_mod[,c('compare', 'compare_num', 'ratio'):=NULL]
rat_table_mod[, vaccination:='No vacc']

all_summary = rbind(oi_summary, rat_table_mod)

all_summary[,Study := suscrule]

vacc_imp  =
  
  ggplot(all_summary[suscrule != '1. Equal']) + 
  geom_pointrange(aes(x =vaccination, y=rat50, ymin=rat05, ymax=rat95, color=Study)) + 
  xlab('Vaccination scenario') + 
  ylab('Change in R relative to Lockdown 3 \n upon reopening of schools')+ 
  scale_y_continuous(
    
    # Features of the first axis
    name = "Relative change in R from Lockdown 3",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*0.8, name="Absolute change from R = 0.8")
  )

ggsave('outputs/plots/Figure1a.png', vacc_imp, width = 15, height = 10)

all_summary_slim = all_summary[suscrule != '1. Equal',c('suscrule', 'vaccination', 'overallimpact')]

all_summary_fat = dcast(all_summary_slim, suscrule ~ vaccination, value.var = 'overallimpact')

