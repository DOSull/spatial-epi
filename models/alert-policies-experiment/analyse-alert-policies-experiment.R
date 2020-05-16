library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/code/spatial-epi/results/alert-policies-experiment")

headers <- list.files(pattern='.header')
casefiles <- list.files(pattern='.cases')

data.all <- data.frame()

for (h in headers) {
  hname <- as.character(h)
  header <- read.csv(hname, skip=2, col.names=c('name', 'value'), header=FALSE)
  csvname <- as.character(header[header$name=='log.file.name', 2])
  
  # should add a function to augment the data with selected header
  # (i.e. global run level) variables
  new_data <- read.csv(csvname) %>% 
    mutate(alert.policy = header[header$name=='alert.policy', 2],
           seed = as.integer(as.character(header[header$name=='seed', 2])))
  data.all <- data.all %>%
    bind_rows(new_data)
}

data.all <- data.all %>%
  filter(ticks <= 300)

data.totals <- data.all %>%
  select(ticks, who, seed, alert.policy, new.cases, alert.level) %>%
  dplyr::group_by(ticks, seed, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.cases), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  dplyr::group_by(ticks, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.cases), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  arrange(alert.policy)


ggplot(data.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=new.cases_min, ymax=new.cases_max), alpha=0.35) +
  geom_line(aes(y=new.cases_mean), lwd=0.2) +
  facet_wrap(~ alert.policy, nrow=2, labeller=(.cols=label_both)) +
  xlab('Days') +
  ylab('Cases per day')

ggsave('cases-time-series-by-alert-policy.png')


all.tick.alert.policy.combos <- crossing(ticks=0:max(data.all$ticks),
                                         alert.level=1:4,
                                         alert.policy=unique(data.all$alert.policy),
                                         seed=1:30)


data.alert.level.totals <- data.all %>%
  select(ticks, seed, alert.level, alert.policy, pop.0) %>%
  dplyr::group_by(ticks, seed, alert.level, alert.policy) %>%
  summarise(total.pop = sum(pop.0), 
            max.ticks = max(ticks)) %>%
  merge(all.tick.alert.policy.combos, all.y=TRUE) %>%
  replace(is.na(.), 0) %>%
  filter(ticks <= max.ticks) %>%
  arrange(ticks)



ggplot(data.alert.level.totals, aes(x=ticks)) +
  geom_step(aes(y=total.pop, group=seed), lwd=0.1, alpha=0.75) +
  facet_grid(alert.level ~ alert.policy, 
             labeller=labeller(.rows=label_both, .cols=label_both)) +
  ylab('Population') +
  xlab('Days')

ggsave('population-in-alert-levels-by-alert-policy.png')
