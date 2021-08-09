ggplot() + 
  geom_line(data=england_rs[between(date, lubridate::ymd('20200901'), lubridate::ymd('20210118'))], aes(x=date, y=mean)) + 
  geom_ribbon(data=england_rs[between(date, lubridate::ymd('20200901'), lubridate::ymd('20210118'))], aes(x=date, ymin=lower_90, ymax=upper_90), fill='red', alpha=0.5)+
  annotate("rect", xmin=lubridate::ymd(20201105), xmax=lubridate::ymd(20201202), ymin=0.4, ymax=2,alpha = .2)+
  annotate("rect", xmin=lubridate::ymd(20210105), xmax=lubridate::ymd(20210118), ymin=0.4, ymax=2,alpha = .2)+
  geom_hline(yintercept = 1.0, linetype='dashed')+
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_date(breaks='month', date_labels = "%b")+
  ylab('R estimate') 


consensus_Rs = fread('../comix/data/consensus_rs.csv')
consensus_Rs[,date:=as.Date(consensus_Rs$date, format = "%d-%b-%y")]
consensus_Rs[,date:= date-10]
ggplot() +
  geom_ribbon(data=consensus_Rs[between(date, lubridate::ymd('20200901'), lubridate::ymd('20210131'))], aes(x=date, ymin=`Lower bound`, ymax=`Upper bound`), fill='red', alpha=0.5)+
  annotate("rect", xmin=lubridate::ymd(20201105), xmax=lubridate::ymd(20201202), ymin=0.4, ymax=2,alpha = .2)+
  annotate("rect", xmin=lubridate::ymd(20210105), xmax=lubridate::ymd(20210213), ymin=0.4, ymax=2,alpha = .2)+
  geom_hline(yintercept = 1.0, linetype='dashed')+
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_date(breaks='month', date_labels = "%b")+
  ylab('R estimate') 


UK_cases = covidregionaldata::get_regional_data(country = "UK")
UK_cases = data.table(UK_cases)
eng_cases = UK_cases[!(region %in% c('Scotland', 'Wales')), sum(cases_new) , by=date]
eng_cases = unique(eng_cases)
names(eng_cases) = c('date', 'cases_new')

ggplot() +
  geom_col(data=eng_cases[between(date, lubridate::ymd('20200901'), lubridate::ymd('20210131'))], aes(x=date, y=cases_new), fill='red', alpha=0.5)+
  annotate("rect", xmin=lubridate::ymd(20201105), xmax=lubridate::ymd(20201202), ymin=0.4, ymax=1.5e5,alpha = .2)+
  annotate("rect", xmin=lubridate::ymd(20210105), xmax=lubridate::ymd(20210131), ymin=0.4, ymax=1.5e5,alpha = .2)+
  scale_x_date(breaks='month', date_labels = "%b")+
  ylab('cases') 

