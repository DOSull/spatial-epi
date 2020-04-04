library(sf)

setwd("~/Documents/code/spatial-epi/geo")

dir()

dhbs <- st_read("nz-district-health-boards-2012.shp")
dhbs.c <- st_centroid(dhbs)

xy <- st_coordinates(dhbs.c)

dhbs.c <- dhbs.c %>% 
  mutate(x = xy[,1],
         y = xy[,2]) %>%
  mutate_each(list(rel=scale), x, y)

pops <- read.csv('dhb-pops.csv')

dhbs.c <- dhbs.c %>% merge(pops, by.x='NAME', by.y='DHB.locality')
dhbs.t <- dhbs.c %>%
  select(NAME, x_rel, y_rel, pop) %>%
  st_drop_geometry()
write.table(dhbs.t, 'dhbs-nl.csv', row.names=FALSE, sep=' ')
