library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/code/spatial-epi/results/random-locales-num-locales")

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
    mutate(num.locales = header[header$name=='num.locales', 2],
           seed = as.integer(as.character(header[header$name=='seed', 2])))
  data.all <- data.all %>%
    bind_rows(new_data)
}

# data.all <- data.all %>%
#   filter(ticks <= 300)

data.totals <- data.all %>%
  select(ticks, who, seed, num.locales, new.cases, alert.level) %>%
  dplyr::group_by(ticks, seed, num.locales) %>%
  summarise_at(vars(-group_cols(), new.cases), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  dplyr::group_by(ticks, num.locales) %>%
  summarise_at(vars(-group_cols(), new.cases), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  mutate(num.locales = as.integer(as.character(num.locales))) %>%
  arrange(num.locales)


ggplot(data.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=new.cases_min, ymax=new.cases_max), alpha=0.35) +
  geom_line(aes(y=new.cases_mean), lwd=0.2) +
  facet_wrap(~ num.locales, nrow=1, labeller=(.cols=label_both)) +
  xlab('Days') +
  ylab('New cases per day')

ggsave('cases-time-series-by-num-locales.png')


all.tick.num.locales.combos <- crossing(ticks=0:max(data.all$ticks),
                                         alert.level=1:4,
                                         alert.policy=unique(data.all$num.locales),
                                         seed=1:30)


data.num.locales.totals <- data.all %>%
  select(ticks, seed, num.locales, alert.level, pop.0) %>%
  dplyr::group_by(ticks, seed, num.locales, alert.level) %>%
  summarise(total.pop = sum(pop.0), 
            max.ticks = max(ticks)) %>%
  merge(all.tick.num.locales.combos, all.y=TRUE) %>%
#  replace(is.na(.), 0) %>%
  filter(ticks <= max.ticks) %>%
  ungroup() %>%
  mutate(num.locales = as.integer(as.character(num.locales))) %>%
  arrange(ticks, num.locales)



ggplot(data.num.locales.totals, aes(x=ticks)) +
  geom_step(aes(y=total.pop, group=seed), lwd=0.1, alpha=0.75) +
  facet_grid(alert.level ~ num.locales, 
             labeller=labeller(.rows=label_both, .cols=label_both)) +
  ylab('Population') +
  xlab('Days')

ggsave('population-in-alert-levels-by-num-locales.png')
