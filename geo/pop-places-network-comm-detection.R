library(sf)
library(tmap)
library(dplyr)
library(dbscan)
library(igraph)
library(ggplot2)

setwd("~/Documents/code/spatial-epi/geo")

roads <- st_read("data/SI-roads.shp")

pop_places <- st_read("data/south-island-populated-places-points.shp")

xy <- st_coordinates(pop_places)
pop_places.m <- pop_places %>%
  mutate(x = xy[, 1],
         y = xy[, 2],
         pop = as.integer(as.character(pop)))

pop_places.id <- pop_places.m %>%
  st_drop_geometry() %>%
  select(gid, name, pop, x, y)

distances <- read.csv("data/SI-road-network-distances.csv", sep=';')

# --------------------------------------------------------
# Building a graph
distances.m <- distances %>%
  na.omit() %>%
  # remove self loops
  filter(origin_id != destination_id) %>%
  # sort by distance from origins
  group_by(origin_id) %>%                                
  arrange(origin_id, network_cost) %>%
  mutate(rank = row_number(), .by_group=TRUE) %>%
  # and retain only the nearest k
  filter(rank < 6) %>%
  ungroup()

edges <- distances.m %>% 
  rename(from=origin_id, to=destination_id, weight=network_cost)
# vertices <- 

G <- graph_from_data_frame(edges, directed=FALSE, vertices=pop_places.id)

plot(G, layout=xy, vertex.size=3, vertex.color='white',
     vertex.label.cex=.5, vertex.label.color='black', vertex.size='pop',
     edge.color='grey')


tmap_mode('view')
tm_shape(pop_places.m) + 
  tm_dots()

