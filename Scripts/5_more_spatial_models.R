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

# 1. Obtain  data ####
pcod_df  <-  sdmTMB::pcod
pcod_sf <- st_as_sf(pcod_df, coords = c("lon","lat"), crs = 4326) %>% 
  st_transform(pcod_sf,
               crs = "+proj=utm +zone=9 +datum=WGS84 +no_defs +type=crs +units=km" ) %>% 
  select(present)


qcs_grid <- sdmTMB::qcs_grid
covars <- rast(qcs_grid, type = "xyz")
crs(covars) <- crs(pcod_sf)

plot(covars)                    
covars_new <- subset(covars, 1)
covars_new$slope <- terrain(covars_new$depth, v = "slope")
covars_new$tri <- terrain(covars_new$depth, v = "TRI", neighbors = 8)

plot(covars_new)
layerCor(covars_new, "cor")

covars_scaled <- scale(covars_new)

plot(covars_scaled)

ggplot() + 
  geom_spatraster(data = covars_scaled, aes(fill = depth)) + 
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = pcod_sf, aes(col = factor(present))) + 
  scale_color_manual(values = c("black","orange"),
                     labels= c("Absence","Presence")) + 
  theme_bw()

# 2. Build mesh ####

mesh <- fm_mesh_2d(loc = pcod_sf, # instead of boundary, we give points
                   cutoff = 2,
                   max.edge = c(7,20), # The largest allowed triangle edge length.
                   offset = c(5,50), 
                   crs = st_crs(pcod_sf)) # The automatic extension distance

ggplot() + 
  gg(mesh) +
  # geom_sf(data= pcod_sf, aes(color = factor(present)), size = 0.5) +
  xlab("") + ylab("")

# 3. Build SPDE ####
st_bbox(pcod_sf)

spde_model <- inla.spde2.pcmatern(mesh,
                                  prior.range = c(100, 0.5),
                                  prior.sigma = c(1, 0.1))


# 4. Components ####

cmp <- ~ Intercept(1) + 
  Eff.depth(covars_scaled$depth, model = "linear") + 
  Eff.slope(covars_scaled$slope, model = "linear") +
  Eff.tri(covars_scaled$tri, model = "linear") +
  Eff.space(geometry, model = spde_model)

# 5. Run model ####

formula <-  present ~ .

lik1 <- bru_obs(formula = formula,
                data = pcod_sf,
                family = "binomial")

m1 <- bru(cmp, 
          lik1)

summary(m1)

plot(spde.posterior(m1, "Eff.space", what = "range")) +
plot(spde.posterior(m1, "Eff.space", what = "log.variance"))

## 5.1. Predict ####

msk <- st_as_sf(st_convex_hull(st_union(pcod_sf)))

new_data <- fm_pixels(mesh, 
                      mask = msk, 
                      dims = c(120, 100)) %>% 
  st_set_crs(st_crs(pcod_sf))

new_rs <- rast(st_coordinates(new_data), type = "xyz")

pred1 <- predict(m1, 
                 new_data,
                 mask = msk, 
                 ~ plogis(Intercept + Eff.depth + Eff.depth + Eff.tri + Eff.space))


pred1_rs <- rast(cbind(st_coordinates(pred1), pred1$q0.5), type = "xyz")

ggplot() + 
  geom_spatraster(data = pred1_rs, aes(fill = X)) +
  scale_fill_viridis_c(na.value = NA) + 
  # gg(mesh) + 
  geom_sf(data = pcod_sf, aes(col = factor(present)), size = 0.5) + 
  scale_color_manual(values = c("black","orange"),
                     labels= c("Absence","Presence")) + 
  theme_bw()

ggplot() + 
  geom_sf(data = pred1, aes(col = q0.5)) +
  scale_colour_viridis_c() +
  ggtitle("Posterior median") +
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = pred1, aes(col = q0.975-q0.025)) +
  scale_colour_viridis_c() +
  ggtitle("Posterior 95% CI") + 
  theme_bw()

## 5.2 Covariate effect

depth.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    depth_new = seq(min(covars_new$depth[], na.rm = T), 
                    max(covars_new$depth[], na.rm = T), 
                    length.out = 100)),
  formula = ~ Eff.depth_eval(depth_new)) 

slope.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    slope_new = seq(min(covars_new$slope[], na.rm = T), 
                    max(covars_new$slope[], na.rm = T), 
                    length.out = 100)),
  formula = ~ Eff.slope_eval(slope_new)) 

tri.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    tri_new = seq(min(covars_new$tri[], na.rm = T), 
                  max(covars_new$tri[], na.rm = T), 
                  length.out = 100)),
  formula = ~ Eff.tri_eval(tri_new)) 

ggplot(depth.pred) +
  geom_line(aes(depth_new, q0.5)) +
  geom_ribbon(aes(depth_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() +

ggplot(slope.pred) +
  geom_line(aes(slope_new, q0.5)) +
  geom_ribbon(aes(slope_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 

ggplot(tri.pred) +
  geom_line(aes(tri_new, q0.5)) +
  geom_ribbon(aes(tri_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()
