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

# 1. Obtaining data ####

data("gorillas_sf")

nests <- gorillas_sf$nests 

counts <- gorillas_sf$plotsample$counts %>% 
  select(-exposure)

mesh <- gorillas_sf$mesh

covars <- rast("Data/covars.grd")
covars_s <- scale(covars)
plot(covars_s)

# 2. Model ####
matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(2, 0.1),
                              prior.sigma = c(0.1, 0.1))
cmp <- ~ -1 +
  spatial_field(geometry, model = matern) +
  inter_cp(1) + inter_pois(1) + scale(1) + 
  Eff.elev(covars_s$elev, model = "linear") + 
  Eff.slope(covars_s$slopeangle, model = "linear")

lik1 <- bru_obs(formula = geometry ~ -1 + inter_cp + spatial_field + 
                  Eff.elev + Eff.slope, 
                family = "cp",
                data = nests,
                domain =  list(geometry = mesh))

lik2 <- bru_obs(formula = count ~ -1 + inter_pois + spatial_field*scale + 
                  Eff.elev + Eff.slope, 
                family = "poisson",
                data = counts,
                domain =  list(geometry = mesh))

m1 <- bru(cmp, lik1, lik2, 
          options = list(control.inla = list(int.strategy = "eb")))

summary(m1)

# 3. Predict ####

newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = gorillas_sf$boundary,
                   format = "sf")

pred <- predict(
  m1, 
  newdf, 
  formula = ~ exp(inter_cp + inter_pois + spatial_field*scale + 
                    Eff.elev + Eff.slope))


ggplot() + 
  gg(data = pred, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  geom_sf(data = counts, aes(col = count)) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") +
  theme_bw() 


# 4. Random spatial effect ####
rdm <- predict(
  m1, 
  newdf, 
  formula = ~ spatial_field)

ggplot() + 
  gg(data = rdm, aes(fill = q0.5), geom = "tile") +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm$q0.5)) * c(-1, 1)) + 
  theme_bw()  
  


bru_names(m1)

plot(m1, "scale")
