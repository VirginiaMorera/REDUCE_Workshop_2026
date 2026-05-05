# 0. Housekeeping ####
rm(list = ls())
load("Environments/4_advanced_spatial_models.RData")

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

# 1. Load and explore data ####
data("gorillas_sf")
names(gorillas_sf)

ggplot() + 
  geom_sf(data = gorillas_sf$boundary, fill = "lightgray") +
  geom_sf(data = gorillas_sf$nests, alpha = 0.2) + 
  geom_sf(data = gorillas_sf$plotsample$counts, aes(col = count), size = 3) +
  scale_colour_viridis_c(option = "B") + 
  theme_bw()

plotsamples <- gorillas_sf$plotsample$counts

# 2. Covars ####
covars <- rast("Data/covars.grd")
covars_s <- scale(covars)
plot(covars_s)

ggplot() + 
  geom_spatraster(data = covars_s, aes(fill = elev)) +
  scale_fill_viridis_c() +
  geom_sf(data = plotsamples, aes(col = count, size = exposure)) +
  scale_colour_viridis_c(option = "B") + 
  theme_bw()

# 3. Mesh ####

ggplot() + 
  gg(gorillas_sf$mesh) + 
  geom_sf(data = gorillas_sf$boundary, col = "red", fill = NA) + 
  geom_sf(data = plotsamples, aes(col = count), size = 3) +
  geom_sf(data = gorillas_sf$nests, alpha = 0.5, size = 1) + 
  scale_colour_viridis_c(option = "B") +
  theme_bw()

outer_boundary <- st_buffer(gorillas_sf$boundary, 2)
inner_boundary <- st_buffer(gorillas_sf$boundary, 0.15)

mesh_better <- fm_mesh_2d_inla(
  boundary = list(inner_boundary, outer_boundary),
  max.edge = c(0.1,0.5),  # this controls the size of the triangles
  crs = st_crs(gorillas_sf$nests)) 

mesh <- gorillas_sf$mesh

ggplot() + 
  gg(mesh_better) + 
  geom_sf(data = plotsamples, aes(col = count), size = 3) +
  geom_sf(data = gorillas_sf$nests, alpha = 0.5, size = 1) + 
  scale_colour_viridis_c(option = "B") +
  theme_bw() + 
  
ggplot() + 
  gg(mesh) + 
  geom_sf(data = plotsamples, aes(col = count), size = 3) +
  geom_sf(data = gorillas_sf$nests, alpha = 0.5, size = 1) + 
  scale_colour_viridis_c(option = "B") +
  theme_bw()

# 4. Define SPDE priors ####

st_bbox(gorillas_sf$boundary)

pcmatern <- inla.spde2.pcmatern(mesh_better, 
                                prior.range = c(2, 0.1),                                
                                prior.sigma = c(1, 0.1))

# 5. Run poisson model ####
cmp <- ~  Intercept(1)  +  
  Eff.elevation(covars$elev, model = "linear") + 
  Eff.slope(covars$slopeangle, model = "linear") + 
  Eff.water(covars$waterdist, model = "linear") + 
  Eff.rdm(geometry, model = pcmatern)

lik1 = bru_obs(formula = count ~ Intercept + Eff.elevation + Eff.rdm,
               family = "poisson",
               data = plotsamples, 
               domain =  list(geometry = mesh_better))

m1 <- bru(cmp, 
          lik1, 
          options = list(E = plotsamples$exposure)) 

summary(m1)

## 5.1 Prediction ####
newdf <- fm_pixels(mesh_better,
                   dims = c(100, 100),
                   mask = gorillas_sf$boundary,
                   format = "sf")

pred1 <- predict(m1, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation + Eff.rdm), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() 

## 5.2 Random effect #### 

rdm1 <- predict(m1, newdata = newdf, 
                 ~  Eff.rdm, 
                 n.samples = 100)

ggplot() + 
  gg(data = rdm1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm1$q0.5)) * c(-1, 1)) + 
  theme_bw()

ggplot() + 
  gg(data = rdm1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm1$q0.5)) * c(-1, 1)) + 
  theme_bw() + 

ggplot() + 
  gg(data = pred1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_viridis_c()+ 
  theme_bw() 

## 5.3 Evaluate covariate effect ####
elev.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

## 5.4 Try with non-linear covariate ####

covars_s$elev_group = inla.group(values(covars_s$elev))

cmp2 <- ~  Intercept(1)  +  
  # Eff.elevation(covars_s$elev_group, model = "rw2") + 
  Eff.elevation(covars_s$elev_group, model = "rw2",
                   hyper = list(
                     prec = list(
                       prior = "pc.prec",
                       param = c(3, 0.1)
                     ))) +
  Eff.rdm(geometry, model = pcmatern)

lik2 = bru_obs(formula = count ~ Intercept + Eff.elevation + Eff.rdm,
               family = "poisson",
               data = plotsamples, 
               domain =  list(geometry = mesh_better))

m2 <- bru(cmp2, 
          lik2, 
          options = list(E = plotsamples$exposure)) 

summary(m2)

pred2 <- predict(m2, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation + Eff.rdm), 
                 n.samples = 100)

rdm2 <- predict(m2, newdata = newdf, 
                ~  Eff.rdm, 
                n.samples = 100)

ggplot() + 
  gg(data = pred2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_viridis_c() + 
  scale_colour_viridis_c(option = "B") + 
  theme_bw() + 

ggplot() + 
  gg(data = rdm2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm2$q0.5)) * c(-1, 1)) + 
  theme_bw() 

elev.pred <- predict(
  m2,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev_group[], na.rm = T), 
                        max(covars_s$elev_group[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()


# 6. Run log gaussian cox process model ####

## 6.1 Prep data a bit different 

nests <- gorillas_sf$nests

pcmatern_cp <- inla.spde2.pcmatern(mesh_better, 
                                   prior.range = c(3, 0.1),                                
                                   prior.sigma = c(0.01, 0.01))

cmp_cp <- ~  Intercept(1)  +  
  Eff.elevation(covars_s$elev, model = "linear") +
  Eff.rdm(geometry, model = pcmatern_cp)


## 6.2 Run model ####
lik_cp = bru_obs(formula = geometry ~ Intercept + Eff.elevation + Eff.rdm,
                 family = "cp",
                 data = nests, 
                 # samplers = gorillas_sf$plotsample$plots,
                 domain =  list(geometry = mesh_better))

m3 <- bru(cmp_cp, 
          lik_cp) 

summary(m3)

## 6.3 Obtain prediction ####
pred3 <- predict(m3, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation + Eff.rdm), 
                 n.samples = 100)

rdm3 <- predict(m3, newdata = newdf, 
                 ~ Eff.rdm, 
                 n.samples = 100)


ggplot() + 
  gg(data = pred3, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() +

ggplot() + 
  gg(data = rdm3, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm3$q0.5)) * c(-1, 1)) + 
  theme_bw() 

## 6.4 Evaluate covariate effect ####
elev.pred <- predict(
  m3,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

## 6.5 Try non-linear covariate ####

cmp_cp2 <- ~  Intercept(1)  +
  # Eff.elevation(covars_s$elev_group, model = "rw2") +
  Eff.elevation(covars_s$elev_group, model = "rw2",
                hyper = list(
                  prec = list(
                    prior = "pc.prec",
                    param = c(3, 0.1)
                  ))) +
  Eff.rdm(geometry, model = pcmatern_cp)

m4 <- bru(cmp_cp2, 
          lik_cp) 

summary(m4)

pred4 <- predict(m4, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation + Eff.rdm), 
                 n.samples = 100)

rdm4 <- predict(m4, newdata = newdf, 
                ~ Eff.rdm, 
                n.samples = 100)

ggplot() + 
  gg(data = pred4, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() +
  
ggplot() + 
  gg(data = rdm4, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm3$q0.5)) * c(-1, 1)) + 
  theme_bw() 

elev.pred <- predict(
  m4,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()
