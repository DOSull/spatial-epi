library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/code/spatial-epi/results/NZ-static-local-scripted-SEIR")

headers <- list.files(pattern='.header')
# casefiles <- list.files(pattern='.cases')

data.all <- data.frame()

for (h in headers) {
  hname <- as.character(h)
  header <- read.csv(hname, skip=2, col.names=c('name', 'value'), header=FALSE)
  csvname <- as.character(header[header$name=='log.file.name', 2])
  
  # should add a function to augment the data with selected header
  # (i.e. global run level) variables
  new_data <- read.csv(csvname) %>% 
    mutate(alert.policy = header[header$name=='alert.policy', 2],
           setup.method = header[header$name=='setup.method', 2],
           start.time = as.integer(as.character(header[header$name=='start.time', 2])),
           seed = as.integer(as.character(header[header$name=='seed', 2])))
  data.all <- data.all %>%
    bind_rows(new_data)
}

dhb.data <- data.all %>%
  filter(setup.method == 'NZ DHBs random cases') %>%
  mutate(t = ticks - start.time)

dhb.totals <- dhb.data %>%
  select(t, who, seed, alert.policy, new.infected, alert.level) %>%
  dplyr::group_by(t, seed, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.infected), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  dplyr::group_by(t, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.infected), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  arrange(alert.policy)


ggplot(dhb.totals, aes(x=t)) +
  geom_ribbon(aes(ymin=new.infected_min, ymax=new.infected_max), alpha=0.35) +
  geom_line(aes(y=new.infected_mean), lwd=0.2) +
  facet_wrap(~ alert.policy, nrow=2, labeller=(.cols=label_both)) +
  xlab('Days') +
  ylab('New symptomatic cases per day') + 
  ggtitle('SEIR model, NZ DHBs, 2000 random cases')

ggsave('seir-nz-dhbs-infected-time-series-by-alert-policy.png')

all.t.alert.policy.combos <- crossing(t=0:max(dhb.data$t),
                                         alert.level=1:4,
                                         alert.policy=unique(dhb.data$alert.policy),
                                         seed=1:30)


dhb.alert.level.totals <- dhb.data %>%
  select(t, seed, alert.level, alert.policy, pop.0) %>%
  dplyr::group_by(t, seed, alert.level, alert.policy) %>%
  summarise(total.pop = sum(pop.0), 
            max.t = max(t)) %>%
  merge(all.t.alert.policy.combos, all.y=TRUE) %>%
  replace(is.na(.), 0) %>%
  filter(t <= max.t) %>%
  arrange(t)

ggplot(dhb.alert.level.totals, aes(x=t)) +
  geom_step(aes(y=total.pop, group=seed), lwd=0.1, alpha=0.75) +
  facet_grid(alert.level ~ alert.policy, 
             labeller=labeller(.rows=label_both, .cols=label_both)) +
  ylab('Population') +
  xlab('Days') +
  ggtitle('SEIR model, NZ DHBs, 2000 random cases')

ggsave('seir-nz-dhbs-population-in-alert-levels-by-alert-policy.png')



ta.data <- data.all %>%
  filter(setup.method == 'NZ TAs random cases') %>%
  mutate(t = ticks - start.time)

ta.totals <- ta.data %>%
  select(t, who, seed, alert.policy, new.infected, alert.level) %>%
  dplyr::group_by(t, seed, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.infected), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  dplyr::group_by(t, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.infected), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  arrange(alert.policy)


ggplot(ta.totals, aes(x=t)) +
  geom_ribbon(aes(ymin=new.infected_min, ymax=new.infected_max), alpha=0.35) +
  geom_line(aes(y=new.infected_mean), lwd=0.2) +
  facet_wrap(~ alert.policy, nrow=2, labeller=(.cols=label_both)) +
  xlab('Days') +
  ylab('New symptomatic cases per day') + 
  ggtitle('SEIR model, NZ TAs, 2000 random cases')

ggsave('seir-nz-tas-infected-time-series-by-alert-policy.png')

all.t.alert.policy.combos <- crossing(t=0:max(ta.data$t),
                                      alert.level=1:4,
                                      alert.policy=unique(ta.data$alert.policy),
                                      seed=1:30)


ta.alert.level.totals <- ta.data %>%
  select(t, seed, alert.level, alert.policy, pop.0) %>%
  dplyr::group_by(t, seed, alert.level, alert.policy) %>%
  summarise(total.pop = sum(pop.0), 
            max.t = max(t)) %>%
  merge(all.t.alert.policy.combos, all.y=TRUE) %>%
  replace(is.na(.), 0) %>%
  filter(t <= max.t) %>%
  arrange(t)

ggplot(ta.alert.level.totals, aes(x=t)) +
  geom_step(aes(y=total.pop, group=seed), lwd=0.1, alpha=0.75) +
  facet_grid(alert.level ~ alert.policy, 
             labeller=labeller(.rows=label_both, .cols=label_both)) +
  ylab('Population') +
  xlab('Days') +
  ggtitle('SEIR model, NZ TAs, 2000 random cases')

ggsave('seir-nz-tas-population-in-alert-levels-by-alert-policy.png')



eradication.times <- data.all %>%
  mutate(t = ticks - start.time) %>%
  select(t, who, seed, alert.policy, setup.method) %>%
  dplyr::group_by(seed, alert.policy, setup.method) %>%
  summarise_at(vars(-group_cols(), t), max) %>%
  ungroup() %>%
  arrange(alert.policy, setup.method)

ggplot(eradication.times, aes(x=alert.policy)) +
  geom_boxplot(aes(y=t)) +
  facet_wrap( ~ setup.method, nrow=1) +
  xlab('Alert policy') + 
  ylab('Time to eradication') +
  ggtitle('SEIR model, 2000 random cases')

ggsave('seir-nz-eradication-times.png')

