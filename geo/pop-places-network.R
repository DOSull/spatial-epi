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
  scale_color_viridis_d(option='D') +
  scale_alpha_continuous()
