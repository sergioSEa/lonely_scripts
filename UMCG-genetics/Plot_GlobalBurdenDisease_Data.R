#Introduction Plots
library(tidyverse)
#GBD Results tool:
#Use the following to cite data included in this download:
#  Global Burden of Disease Collaborative Network.
#Global Burden of Disease Study 2019 (GBD 2019) Results.
#Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2020.
#Available from https://vizhub.healthdata.org/gbd-results/.


read_csv("~/Downloads/IHME-GBD_2019_DATA-bce1875a-1/IHME-GBD_2019_DATA-bce1875a-1.csv") -> Data3
Data3 %>% group_by(year) %>% summarise()
Data3 %>% group_by(cause_name) %>% summarise()

Run_analysis = function(Data ){
  All_stats = tibble()
  for (Disease in unique(Data3$cause_name) ){
    Data %>% filter(cause_name == Disease) %>% filter(metric_name == "Percent") -> data_analysis
    data_analysis %>%  lm(val ~ year, .  ) %>% summary() -> Rs
    Rs$coefficients["year",]  %>% t() %>% as_tibble() %>% mutate(Condition = Disease) %>% rbind(All_stats, .) -> All_stats
    P =Rs$coefficients["year", "Pr(>|t|)"]
    if ( is.na(P) ){ next } 
    print(paste0(Disease," ", P) )
    if (P  > 0.05/13 ){ next }
    data_analysis %>% filter(year  == 1990 | year == 2019 ) %>% print()
    data_analysis %>% ggplot(aes(x=year, y=val )) + geom_point() + theme_bw() + ylab( paste0("Percentage of ", Disease) ) + geom_ribbon( aes(ymin = lower, ymax = upper), fill = "lightgray", alpha = 0.5) %>% print()
  }
  All_stats %>% drop_na() %>% mutate(P_bonf = p.adjust(`Pr(>|t|)`, "bonferroni") ) %>% arrange(P_bonf) %>% return()
  
}
All_stats = Run_analysis( Data3 %>% filter(location_name == "Europe" ) )
All_stats_NL = Run_analysis( Data3 %>% filter(location_name == "Netherlands" ) ) 

c25 <- c("dodgerblue2","green4",  "#E31A1C","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

Data3 %>%  select(cause_name, year, val, upper, lower, location_name) %>% rename(Location = location_name) %>% mutate(Location = ifelse(Location == "Netherlands", "The Netherlands", Location ) ) %>% 
  mutate(cause_name = ifelse(cause_name == "Total burden related to Non-alcoholic fatty liver disease (NAFLD)" , "NAFLD" , cause_name )) %>% 
  mutate(cause_name = ifelse(cause_name == "Endocrine, metabolic, blood, and immune disorders" , "Endocrine, metabolic\nblood and immune" , cause_name ))  %>%
  ggplot(aes(x=year, y= val, col=Location )) + geom_point() + theme_bw() + ylab( paste0("Prevalence of disorder") ) + 
  geom_ribbon( aes(ymin = lower, ymax = upper), fill = "lightgray", alpha = 0.01) + scale_color_manual(values = c25) + geom_line() + facet_wrap(~cause_name, scales = "free") +
  theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20) )  


Data3  %>% filter(cause_name=="Diabetes mellitus type 2") %>% filter(year  == 1990 | year == 2019 ) %>% print()
Data3  %>% filter(cause_name=="Digestive diseases") %>% filter(year  == 1990 | year == 2019 ) %>% print()
Data3  %>% filter(cause_name=="Endocrine, metabolic, blood, and immune disorders") %>% filter(year  == 1990 | year == 2019 ) %>% print()
