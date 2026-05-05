# 0. Housekeeping ####
rm(list = ls())

library(inlabru)
library(INLA)
library(tidyverse)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(patchwork)

bru_options_set(control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       cpo = TRUE))
# 1. Load data ####

badger_data <- readRDS("Data/badger_data.RDS")
carlow <- badger_data$carlow
ca_data_thinned <- badger_data$data
env_vars_scaled <- readRDS("Data/env_vars_ca.RDS")

elevation <- env_vars_scaled$elevation
slope <- env_vars_scaled$slope
forestDist <- env_vars_scaled$forest_distances
topo_wetness <- env_vars_scaled$topographic_wetness_index

# 2. Mesh ####
inner_boundary <- st_buffer(st_simplify(carlow, dTolerance = 1, TRUE),  1)
outer_boundary <- st_buffer(inner_boundary, 15)

mesh <- fm_mesh_2d_inla(
  boundary = list(inner_boundary, outer_boundary),
  max.edge = c(2, 7),  # this controls the size of the triangles
  cutoff = 1, 
  crs = st_crs(ca_data_thinned)) 

ggplot() + 
  gg(mesh) + 
  geom_sf(data = carlow, fill = NA, col = "red") +
  geom_sf(data = ca_data_thinned) + 
  theme_bw()

# 3. Model ####

matern_p <- inla.spde2.pcmatern(mesh,
                                prior.range = c(15, 0.1),
                                prior.sigma = c(1, 0.1))

matern_m <- inla.spde2.pcmatern(mesh,
                                prior.range = c(15, 0.1),
                                prior.sigma = c(1, 0.1))

cmp <- ~ -1 +
  point_field(geometry, model = matern_p) +
  mark_field(geometry, model = matern_m) +
  scale(1) + inter_point(1) + inter_mark(1) +
  Eff.elev_point(elevation, model = "linear") + 
  Eff.slope_point(slope, model = "linear") + 
  Eff.elev_mark(elevation, model = "linear") + 
  Eff.slope_mark(slope, model = "linear") 

lik1 <- bru_obs(formula = geometry ~ -1 + inter_point + point_field + 
                  Eff.elev_point + Eff.slope_point, 
                family = "cp",
                data = ca_data_thinned,
                domain =  list(geometry = mesh))

lik2 <- bru_obs(formula = GROUP_SIZE ~ -1 + inter_mark + point_field*scale + mark_field + 
                  Eff.elev_mark + Eff.slope_mark, 
                family = "poisson",
                data = ca_data_thinned,
                domain =  list(geometry = mesh))

fit <- bru(cmp, lik1, lik2, 
           options = list(control.inla = list(int.strategy = "eb")))
summary(fit)

# 4. Predict ####

newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = carlow,
                   format = "sf")

pred_setts <- predict(
  fit, 
  newdf, 
  formula = ~ exp(inter_point +  point_field + 
                    Eff.elev_point + Eff.slope_point))

pred_badgers <- predict(
  fit, 
  newdf, 
  formula = ~ exp(inter_mark +  point_field*scale + 
                    Eff.elev_mark + Eff.slope_mark))


ggplot() + 
  gg(data = pred_setts, aes(fill = q0.5), geom = "tile") +
  # geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() + 
  
  ggplot() + 
  gg(data = pred_badgers, aes(fill = q0.5), geom = "tile") +
  # geom_sf(data = nests, aes(col = mark), size = 1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() 
