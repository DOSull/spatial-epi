library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/code/spatial-epi/results/NZ-historical")

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
           setup.method = header[header$name=='setup.method', 2],
           start.time = as.integer(as.character(header[header$name=='start.time', 2])),
           seed = as.integer(as.character(header[header$name=='seed', 2])))
  data.all <- data.all %>%
    bind_rows(new_data)
}

dhb.data <- data.all %>%
#  filter(setup.method == 'NZ DHBs random cases') %>%
  mutate(t = ticks - start.time)

dhb.totals <- dhb.data %>%
  select(t, who, seed, alert.policy, new.cases, alert.level) %>%
  dplyr::group_by(t, seed, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.cases), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  dplyr::group_by(t, alert.policy) %>%
  summarise_at(vars(-group_cols(), new.cases), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  arrange(alert.policy)


ggplot(dhb.totals, aes(x=t)) +
  geom_ribbon(aes(ymin=new.cases_min, ymax=new.cases_max), alpha=0.35) +
  geom_line(aes(y=new.cases_mean), lwd=0.2) +
  facet_wrap(~ alert.policy, nrow=1, labeller=(.cols=label_both)) +
  xlab('Days') +
  ylab('New symptomatic cases per day') + 
  ggtitle('Branching model, NZ DHBs, 300 random cases')

ggsave('branching-nz-dhbs-new-cases-time-series-by-alert-policy.png')




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
  ggtitle('Branching model, NZ DHBs, 300 random cases')

ggsave('branching-nz-dhbs-population-in-alert-levels-by-alert-policy.png')






