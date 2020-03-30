devtools::install_github("terminological/uk-covid-datatools", force = TRUE)

library(ukcovidtools)
library(ggplot2)
library(EpiEstim)
library(tidyverse)

#UK = ukcovidtools::getUKCovidTimeseries()
#cleanRegional = UK$tidyUKRegional %>% ukcovidtools::normaliseAndCleanse()

# A list of 7 data sets of similar nature and formats
#glimpse(UK$UKregional)
#glimpse(UK$englandNHS)
#glimpse(UK$englandUnitAuth)
#glimpse(UK$englandUnitAuth2NHSregion)
#glimpse(UK$tidyUKRegional)
#glimpse(UK$tidyEnglandNHS)
#glimpse(UK$tidyEnglandUnitAuth)


#siConfig = EpiEstim::make_config(list(
#  mean_si = 4.7, 
#  std_si = 2.9
#))

#regionalRt = cleanRegional %>% ukcovidtools::tidyEstimateRt(siConfig)
#ggplot(regionalRt, aes(x=date,y=`Median(R)`,ymin=`Quantile.0.05(R)`,ymax=`Quantile.0.95(R)`,fill=uk_region,colour=uk_region))+
#  geom_ribbon(alpha=0.2, colour=NA)+geom_line()+geom_hline(yintercept = 1, colour="red")+expand_limits(y=0)


#### Lockdown impact ####
library(rgdal); library(ggplot2); library(ggspatial); library(rgeos); library(maptools); library(lubridate)
library(patchwork); library(sp); library(ggrepel); library(gganimate); library(gifski); library(transformr)
#ggplot2::theme_set(standardPrintOutput::defaultFigureLayout())

serialIntervals = tibble(
  mean_si_estimate = c(3.96, 6.3, 4.22, 4.56, 3.95, 5.21, 4.7, 7.5,6.6),
  mean_si_estimate_low_ci = c(3.53, 5.2, 3.43, 2.69,-4.47, -3.35, 3.7, 5.3, 0.7),
  mean_si_estimate_high_ci = c(4.39, 7.6, 5.01, 6.42, 12.51,13.94, 6.0, 19.0, 19.0),
  std_si_estimate = c(4.75,4.2, 0.4, 0.95, 4.24, 4.32, 2.3, 3.4, NA),
  std_si_estimate_low_ci = c(4.46, 3.1, NA, NA, 4.03, 4.06, 1.6, NA, NA),
  std_si_estimate_high_ci = c(5.07, 5.3, NA, NA, 4.95, 5.58, 3.5, NA, NA),
  sample_size = c(468,48,135,93,45,54,28,16,90),
  population = c("China", "Shenzhen","Taijin","Singapore","Taijin","Singapore", "SE Asia", "Wuhan","Italy"),
  source = c(
    "Zhanwei Du et al. Serial Interval of COVID-19 among Publicly Reported Confirmed Cases. Emerging Infectious Disease journal 26, (2020)",
    "Bi, Q. et al. Epidemiology and Transmission of COVID-19 in Shenzhen China: Analysis of 391 cases and 1,286 of their close contacts. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.03.20028423",
    "Tindale, L. et al. Transmission interval estimates suggest pre-symptomatic spread of COVID-19. Epidemiology (2020) doi:10.1101/2020.03.03.20029983",
    "Tindale, L. et al. Transmission interval estimates suggest pre-symptomatic spread of COVID-19. Epidemiology (2020) doi:10.1101/2020.03.03.20029983",
    "Ganyani, T. et al. Estimating the generation interval for COVID-19 based on symptom onset data. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.05.20031815",
    "Ganyani, T. et al. Estimating the generation interval for COVID-19 based on symptom onset data. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.05.20031815",
    "Nishiura, H., Linton, N. M. & Akhmetzhanov, A. R. Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. (2020) doi:10.1016/j.ijid.2020.02.060",
    "Li, Q. et al. Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus-Infected Pneumonia. N. Engl. J. Med. (2020) doi:10.1056/NEJMoa2001316",
    "Cereda, D. et al. The early phase of the COVID-19 outbreak in Lombardy, Italy. arXiv [q-bio.PE] (2020)")
)
unk=function(x) ifelse(is.na(x),"unk",x)
table = serialIntervals %>% mutate(
  `Mean SI\n(95% CrI) days`=paste0(mean_si_estimate,"\n(",unk(mean_si_estimate_low_ci),"-",
                                   unk(mean_si_estimate_high_ci),")"),
  `Std SI\n(95% CrI) days`=paste0(unk(std_si_estimate),"\n(",unk(std_si_estimate_low_ci),"-",unk(std_si_estimate_high_ci),")")
) %>% select(-contains("estimate")) %>% select(
  `Reference`=source,
  `Mean SI\n(95% CrI) days`,
  `Std SI\n(95% CrI) days`,
  `N`=sample_size,
  `Population`=population
)
#table %>% group_by(`Reference`) %>% standardPrintOutput::saveTable("~/Dropbox/covid19/lockdown-impact/serialIntervals", defaultFontSize = 8, colWidths=c(4.5,2,2,0.5,1))

wtSIs = serialIntervals %>% summarise(
  mean_si = weighted.mean(mean_si_estimate,sample_size,na.rm = TRUE),
  min_mean_si = weighted.mean(mean_si_estimate_low_ci,sample_size,na.rm = TRUE),
  max_mean_si = weighted.mean(mean_si_estimate_high_ci,sample_size,na.rm = TRUE),
  std_si  = weighted.mean(ifelse(is.na(std_si_estimate_low_ci),NA,1)*std_si_estimate,sample_size,na.rm = TRUE),
  min_std_si  = weighted.mean(std_si_estimate_low_ci,sample_size,na.rm = TRUE),
  max_std_si  = weighted.mean(std_si_estimate_high_ci,sample_size,na.rm = TRUE)
  #total = sum(sample_size)
) %>% mutate(
  std_mean_si = (max_mean_si - min_mean_si) / 3.92, # TODO: fit gamma
  std_std_si = (max_std_si - min_std_si) / 3.92
)
tdp = function(x,y,z) sprintf("%1.2f (%1.2f - %1.2f)", x ,y, z)
tdp(wtSIs$mean_si, wtSIs$min_mean_si, wtSIs$max_mean_si)
tdp(wtSIs$std_si, wtSIs$min_std_si, wtSIs$max_std_si)

keyDates = tibble(
  date = as.Date(c("2020-03-13","2020-03-16","2020-03-19","2020-03-23")), #max(r0shapes$date-1, na.rm=TRUE)),
  impactDate = as.Date(c("2020-03-14","2020-03-21","2020-03-24","2020-03-28")), #max(r0shapes$date-1, na.rm=TRUE)),
  event = c("Inpatient only testing","Social isolation of vulnerable","Travel ban / school closure","Stay at home") #,"Latest")
) %>% mutate(label = paste0(date,": \n",event))

ts = ukcovidtools::getUKCovidTimeseries()
#TODO: update these estimates of SI and use credible interval
# cfg = EpiEstim::make_config(list(
#   mean_si = wtSIs$mean_si,
#   std_si = wtSIs$std_si
# ))
cfg = EpiEstim::make_config(list(
  mean_si = wtSIs$mean_si, 
  std_mean_si = wtSIs$std_mean_si,
  min_mean_si = wtSIs$min_mean_si, 
  max_mean_si = wtSIs$max_mean_si,
  std_si = wtSIs$std_si, 
  std_std_si = wtSIs$std_si,
  min_std_si = wtSIs$min_std_si, 
  max_std_si = wtSIs$max_std_si), method="uncertain_si")

R0regionaltimeseries = ts$tidyUKRegional %>% group_by(uk_region) %>% normaliseAndCleanse() %>% tidyEstimateRt(cfg, window = 5) 
R0regionaltimeseries = R0regionaltimeseries %>% filter(!is.na(`Median(R)`))
ukregionalplot = ggplot(R0regionaltimeseries, aes(x=date, y=`Median(R)`, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`))+ #, colour=uk_region, fill=uk_region))+
  geom_ribbon(alpha=0.2)+geom_line()+geom_hline(yintercept = 1, colour="grey50", linetype="dashed")+facet_wrap(vars(uk_region)) + 
  coord_cartesian(ylim=c(0, 5))+
  geom_vline(aes(xintercept=date,colour=event),data=keyDates, show.legend = FALSE)+
  #geom_vline(aes(xintercept=impactDate,colour=event,linetype="dashed"),data=keyDates, show.legend = FALSE)+
  ggrepel::geom_text_repel(
    aes(x=date, y=Inf, colour=event, label=event),data=keyDates, hjust=0,vjust=1, angle=90, show.legend = FALSE,box.padding=0.05,inherit.aes = FALSE,
    size=(10/ggplot2:::.pt/(96/72)))+
  scale_x_date(date_breaks="1 day", date_labels = "%d-%b")+theme(axis.text.x=element_text(angle = 90, vjust =0.5))
ukregionalplot

R0nhstimeseries = ts$tidyEnglandNHS %>% group_by(england_nhs_region) %>% normaliseAndCleanse() %>% tidyEstimateRt(cfg, window = 5) 
R0nhstimeseries = R0nhstimeseries %>% filter(!is.na(`Median(R)`))
englandnhsplot = ggplot(R0nhstimeseries, aes(x=date, y=`Median(R)`, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`))+ #, colour=uk_region, fill=uk_region))+
  geom_ribbon(alpha=0.2)+geom_line()+
  geom_hline(yintercept = 1, colour="grey50", linetype="dashed")+facet_wrap(vars(england_nhs_region)) + 
  #standardPrintOutput::narrowAndTall()+
  coord_cartesian(ylim=c(0, 5))+
  geom_vline(aes(xintercept=date,colour=event),data=keyDates)+
  #geom_vline(aes(xintercept=impactDate,colour=event,linetype="dashed"),data=keyDates, show.legend = FALSE) #, show.legend = FALSE) #+
  # ggrepel::geom_text_repel(
  #         aes(x=date, y=Inf, colour=event, label=event),data=keyDates, hjust=0,vjust=1, angle=90, show.legend = FALSE,box.padding=0.05,inherit.aes = FALSE,
  #         size=(10/ggplot2:::.pt/(96/72)))
  scale_x_date(date_breaks="1 day", date_labels = "%d-%b")+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
englandnhsplot

R0timeseries = ts$tidyEnglandUnitAuth %>% normaliseAndCleanse() %>% group_by(GSS_CD, GSS_NM) %>% tidyEstimateRt(cfg, window=5)

save(R0timeseries, file = 'R0timeseries.RData')
#write_csv(R0timeseries, "Covid19/Supplementary_Rt_Timeseries_by_Unitary_Authority.csv")

data("UKCovidMaps")

r0shapes = UKCovidMaps$unitaryAuthority %>% 
  left_join(R0timeseries, by=c("ctyua19cd"="GSS_CD")) %>% 
  mutate(ago=difftime(date,lubridate::now(),units="days")) %>% 
  filter(!is.na(date))
r0shapes = r0shapes %>% mutate(`Median(R)` = ifelse(`Median(R)`>10, 9.999,`Median(R)`))
r0shapes_key = r0shapes %>% inner_join(keyDates, by="date")
ukwide = ggplot(r0shapes_key)+
  geom_sf(aes(fill=`Median(R)`),data=r0shapes_key)+
  scale_fill_gradient2(
    low="green",
    mid="white",
    high="red",
    midpoint=0,
    trans="log",
    na.value = "grey80", 
    limits=c(0.1,10), 
    breaks=c(0.1,0.4,1,2.5,10), 
    labels=c("<0.1","0.4","1","2.5",">10"))+
#  standardPrintOutput::narrowAndTall()+
#  standardPrintOutput::mapTheme()+
  facet_wrap(vars(label), nrow = 1)

london = ukwide + coord_sf(crs = 4326,xlim = c(-0.7, 0.5), ylim = c(51.25, 51.75), expand = FALSE)

# layout <- c(
#  patchwork::area(t = 1, l = 1, b = 10, r = 10),
#  patchwork::area(t = 1, l = 1, b = 3, r = 3)
# )
# combined = ukwide / london + plot_layout(guides="collect")
#ukwide %>% standardPrintOutput::saveThirdPageFigure("~/Dropbox/covid19/lockdown-impact/englandMap")
#london %>% standardPrintOutput::saveThirdPageFigure("~/Dropbox/covid19/lockdown-impact/londonMap")
#ukwide
#london

ukwide = ggplot(r0shapes)+
  geom_sf(aes(fill=`Median(R)`), data=r0shapes)+
  scale_fill_gradient2(
    low="green",
    mid="white",
    high="red",
    midpoint=0,
    trans="log",
    na.value = "grey80", 
    limits=c(0.1,10), 
    breaks=c(0.1,0.4,1,2.5,10), 
    labels=c("<0.1","0.4","1","2.5",">10"))
#  standardPrintOutput::narrowAndTall()+
#  standardPrintOutput::mapTheme()
london = ukwide + coord_sf(crs = 4326,xlim = c(-0.7, 0.5), ylim = c(51.25, 51.75), expand = FALSE)

#ukwide
#london

anim = ukwide+gganimate::transition_time(date)
gif = gganimate::animate(anim, renderer=gganimate::gifski_renderer())
gganimate::anim_save("Covid19/Rt_by_unitary_authority_over_time.gif",gif)
anim2 = london+gganimate::transition_time(date)
gif2 = gganimate::animate(anim2, renderer=gganimate::gifski_renderer())
gganimate::anim_save("Covid19/London_Rt_by_unitary_authority_over_time.gif",gif2)







