library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/code/spatial-epi/distributed-seir-results")

header <- read.csv(dir()[1], skip=2, nrows=31, col.names=c('name', 'value'), header=FALSE)
data <- read.csv(dir()[1], skip=33) %>%
  mutate(num.locales = header[header$name=='num-locales', 2],
         seed = header[header$name=='seed', 2])

for (fname in dir()[2:90]) {
  header <- read.csv(fname, skip=2, nrows=31, col.names=c('name', 'value'), header=FALSE)
  add_data <- read.csv(fname, skip=33) %>%
    mutate(num.locales = header[header$name=='num-locales', 2],
           seed = header[header$name=='seed', 2])
  data <- data %>%
    bind_rows(add_data)
}

data.totals <- data %>%
  mutate(num.locales = as.numeric(num.locales)) %>%
  arrange(num.locales) %>%
  group_by(ticks, seed, num.locales, alert.level) %>%
  summarise(pop=sum(pop.0), num=n(), dead=sum(dead), infected=sum(infected), recovered=sum(recovered)) %>%
  ungroup() %>%
  group_by(ticks, alert.level, num.locales) %>%
  mutate(mean.pop = mean(pop), min.pop = min(pop), max.pop = max(pop)) %>%
  ungroup() %>%
  group_by(ticks, num.locales) %>%
  mutate(mean.dead = mean(dead), min.dead = min(dead), max.dead = max(dead),
         mean.infected = mean(infected), min.infected = min(infected), max.infected = max(infected),
         mean.recovered = mean(recovered), min.recovered = min(recovered), max.recovered = max(recovered))


# THIS SEEMS FINE
ggplot(data.m, aes(x=ticks)) +
  geom_ribbon(aes(ymin=min.pop, ymax=max.pop), alpha=0.35, lwd=0) +
  geom_line(aes(y=mean.pop), lwd=0.5) + 
  facet_grid(alert.level ~ num.locales, labeller=labeller(.rows=label_both, .cols=label_both)) 

# THIS IS WEIRD
ggplot(data.m, aes(x=ticks)) +
  geom_ribbon(aes(ymin=min.dead, ymax=max.dead), alpha=0.35, lwd=0) +
  geom_line(aes(y=mean.dead), lwd=0.5) + 
  facet_wrap(~num.locales, nrow=1, labeller=labeller(.cols=label_both))
