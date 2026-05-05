# 0. Housekeeping ####
rm(list = ls())
load("Environments/6_marked_pp.RData")

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
# 1. Read data ####

rc_model <- readRDS("Data/rc_model.RDS")
rc_data <- rc_model$data
rc <- rc_model$rc
covar_rc <- readRDS("Data/covar_rc.RDS.")
plot(covar_rc)

ggplot() + 
  geom_sf(data = rc) + 
  geom_sf(data = rc_data %>% filter(!is.na(Count)), aes(col = Count)) + 
  geom_sf(data = rc_data %>% filter(is.na(Count)), col = "red") +
  theme_bw()

# 2. Create mesh ####

inner_boundary <- st_buffer(st_simplify(rc, dTolerance = 1, TRUE),  1)
outer_boundary <- st_buffer(inner_boundary, 15)

mesh <- fm_mesh_2d_inla(
  boundary = list(inner_boundary, outer_boundary),
  max.edge = c(2, 7),  # this controls the size of the triangles
  cutoff = 1, 
  crs = st_crs(rc_data)) 

ggplot() + 
  gg(mesh) + 
  geom_sf(data = rc, fill = NA, col = "red") +
  geom_sf(data = rc_data) + 
  theme_bw()


# 3. Model ####

matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(30, 0.1),
                              prior.sigma = c(1, 0.1))

cmp <- ~ -1 +
  rdm_field(geometry, model = matern) +
  inter_count(1) + inter_pa(1) + scale(1) + 
  Eff.elev(covar_rc$elevation, model = "linear") + 
  Eff.forest(covar_rc$forest_distances, model = "linear") + 
  Eff.human(covar_rc$human_footprint_index, model = "linear")


lik1 <- bru_obs(formula = Count ~ -1 + inter_count + rdm_field + 
                  Eff.elev + Eff.forest + Eff.human, 
                family = "poisson",
                data = rc_data,
                domain =  list(geometry = mesh))

lik2 <- bru_obs(formula = PA ~ -1 + inter_pa + scale*rdm_field + 
                  Eff.elev + Eff.forest + Eff.human, 
                family = "binomial",
                data = rc_data,
                domain =  list(geometry = mesh))

m1 <- bru(cmp, lik1, lik2, 
          options = list(control.inla = list(int.strategy = "eb")))

summary(m1)

# 4. Predict ####
newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = rc,
                   format = "sf")

pred_prop <- predict(
  m1, 
  newdf, 
  formula = ~ plogis(inter_pa + scale * rdm_field + 
                       Eff.elev + Eff.forest + Eff.human))

pred_count <- predict(
  m1, 
  newdf, 
  formula = ~ exp(inter_count + rdm_field + 
                    Eff.elev + Eff.forest + Eff.human) 
  
)

ggplot() + 
  gg(data = pred_prop, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = rc_data, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw()  + 
  
ggplot() + 
  gg(data = pred_count, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = rc_data, aes(col = Count), size = 1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() 

pred_both <- predict(
  m1, 
  newdf, 
  formula = ~ exp(inter_count + rdm_field + 
                    Eff.elev + Eff.forest + Eff.human)*
    plogis(inter_pa + scale * rdm_field + 
             Eff.elev + Eff.forest + Eff.human)
)

ggplot() + 
  gg(data = pred_both, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = rc_data, aes(col = Count), size = 1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() 

