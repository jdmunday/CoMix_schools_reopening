library(viridis)
library(data.table)
library(ggplot2)
library(ggridges)
library(patchwork)
library(wesanderson)
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
  
  


layout <- "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
"

Fig2 = p1 + p2 + 
  plot_layout(design = layout)

ggsave('outputs/plots/Figure2.pdf', Fig2, width = 15, height = 10)


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
