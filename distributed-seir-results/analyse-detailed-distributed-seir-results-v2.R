library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/code/spatial-epi/distributed-seir-results")

headers <- list.files(pattern='.header')
csvs <- list.files(pattern='.csv')

data <- data.frame()

for (i in 1:length(headers)) {
  header <- read.csv(headers[i], skip=2, col.names=c('name', 'value'), header=FALSE)
  new_data <- read.csv(csvs[i]) %>%
    mutate(num.locales = header[header$name=='num-locales', 2],
           seed = header[header$name=='seed', 2])
  data <- data %>%
    bind_rows(new_data)
}

data.sel <- data %>%
  select(ticks, who, seed, num.locales, pop.0, infected, recovered, dead, alert.level)


data.totals <- data.sel %>%
  group_by(ticks, seed, num.locales) %>%
  summarise_at(vars(-group_cols(), infected, recovered, dead), sum) %>%
  ungroup() %>%
  select(-alert.level) %>%
  group_by(ticks, num.locales) %>%
  summarise_at(vars(-group_cols(), infected, recovered, dead), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  mutate(num.locales = as.numeric(num.locales)) %>%
  arrange(num.locales) %>%
  mutate(num.locales = as.factor(num.locales))
  

p.dead <- ggplot(data.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=dead_min, ymax=dead_max), alpha=0.5)+ 
  geom_line(aes(y=dead_mean)) +
  facet_wrap(~ num.locales, nrow=1, labeller=(.cols=label_both)) +
  ylab('Dead') + 
  xlab('Days')

p.recovered <- ggplot(data.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=recovered_min, ymax=recovered_max), alpha=0.5)+ 
  geom_line(aes(y=recovered_mean)) +
  facet_wrap(~ num.locales, nrow=1, labeller=(.cols=label_both)) +
  ylab('Recovered') + 
  xlab('Days')

p.infected <- ggplot(data.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=infected_min, ymax=infected_max), alpha=0.5)+ 
  geom_line(aes(y=infected_mean)) +
  scale_y_log10() + 
  facet_wrap(~ num.locales, nrow=1, labeller=(.cols=label_both)) +
  ylab('log Infected') + 
  xlab('Days')

ggarrange(p.dead, p.recovered, p.infected, nrow=3)
ggsave('pandemic-time-series-by-num-locales.png')



data.alert.level.totals <- data.sel %>%
  group_by(ticks, seed, num.locales, alert.level) %>%
  summarise_at(vars(-group_cols(), pop.0), sum) %>%
  group_by(ticks, alert.level, num.locales) %>%
  summarise_at(vars(-group_cols(), pop.0), list(~ mean(.), ~ min(.), ~ max(.))) %>%
  ungroup() %>%
  mutate(num.locales = as.numeric(num.locales)) %>%
  arrange(num.locales) %>%
  mutate(num.locales = as.factor(num.locales))

ggplot(data.alert.level.totals, aes(x=ticks)) +
  geom_ribbon(aes(ymin=pop.0_min, ymax=pop.0_max), alpha=0.5) + 
  geom_line(aes(y=pop.0_mean), lwd=0.5) + 
  facet_grid(alert.level ~ num.locales, labeller=labeller(.rows=label_both, .cols=label_both)) +
  ylab('Population') + 
  xlab('Days')

ggsave('population-in-different-alert-levels-by-num-locales.png')


