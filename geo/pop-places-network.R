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

# Make into a matrix
ids <- pop_places.id$gid

o.ids <- match(distances$origin_id, ids)
d.ids <- match(distances$destination_id, ids)

od.ids <- cbind(o.ids, d.ids)

OD.dist <- matrix(0, nrow=length(ids), ncol=length(ids))
OD.dist[od.ids] <- distances$network_cost 
OD.dist[is.na(OD.dist)] <- 1e10
diag(OD.dist) <- 0
distance.matrix <- as.dist(OD.dist)

# DBSCAN clusters on network distances
hdbc <- hdbscan(distance.matrix, minPts=5)
clusters <- data.frame(gid=ids, cluster=hdbc$cluster, prob=hdbc$membership_prob)

pop_places.m <- pop_places.m %>%
  merge(clusters) %>%
  mutate(cluster = as.factor(cluster))

ggplot(roads) +
  geom_sf(col='pink', lwd=0.15) +
  geom_sf(data=pop_places.m, aes(colour=cluster, alpha=prob), cex=3, lwd=0) +
  scale_color_viridis_d(option='D') 
# +
#   geom_sf_text(aes(label=name), cex=0.65)

tm_shape(filter(pop_places.m, prob > 0.25)) +
  tm_dots(col='cluster', palette='Spectral')
# --------------------------------------------------------
# Building a graph
distances.m <- distances %>%
  #  filter(origin_id != destination_id) %>%
  group_by(origin_id) %>%
  arrange(origin_id, network_cost) %>%
  mutate(rank = row_number(), .by_group=TRUE) %>%
  filter(rank < 4) %>%
  ungroup()

G <- graph_from_data_frame(rename(distances.m, from=origin_id, to=destination_id, weight=network_cost), 
                           vertices=pop_places.id, directed=FALSE)

plot(G, layout=xy, vertex.size=3, vertex.color='white',
     vertex.label.cex=.5, vertex.label.color='black', vertex.size='pop',
     edge.color='grey')

