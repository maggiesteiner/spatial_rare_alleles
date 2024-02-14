library(tidyverse)

parents <- read.delim(file="../simulation_checks/nonspatial_parents.txt",header = F,sep = " ") 

colnames(parents) <- c("s","tick","p1","p1_nmut","p2","p2_nmut")

parents %>% filter(tick>70) %>%
  pivot_longer(cols=c(p1,p2),
               names_to = "type",
               values_to = "parent") %>%
  rowwise() %>%
  mutate(num_mut=case_when(type=="p1"~p1_nmut,
                           type=="p2"~p2_nmut)) -> parents_long
# get number of offspring and variance
parents_long %>%
  group_by(tick,parent,carrier=num_mut>0) %>%
  summarise(n_child=n()) %>%
  ungroup() %>%
  group_by(carrier) %>%
  summarise(mn_offspring = mean(n_child),v_offspring = var(n_child)) 
# check age of carriers and noncarriers 
parents_long %>%
  filter(type=="p1") %>% # on p1 because p2 may not be chosen but p1 is guaranteed
  group_by(parent) %>%
  summarise(carrier=num_mut>0,
            tick_appears=min(tick),
            tick_dies=max(tick),
            lifespan=tick_dies-tick_appears + 1) %>% # +1 because chosen for reproduction in the tick after they are born
  filter(tick_appears!=71,tick_dies<100) %>%
  ungroup() %>%
  group_by(carrier) %>%
  summarise(age_at_death=mean(lifespan))
