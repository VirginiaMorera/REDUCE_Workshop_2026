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
covars_new$depth <- covars$depth
covars_new$slope <- terrain(covars_new$depth, v = "slope")
covars_new$tri <- terrain(covars_new$depth, v = "TRI", neighbors = 8)

plot(covars_new)
layerCor(covars_new, "cor")

covars_scaled <- scale(covars_new)

plot(covars_scaled)

ggplot() + 
  geom_spatraster(data = depth, aes(fill = depth_scaled)) + 
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = pcod_sf, aes(col = factor(present))) + 
  scale_color_manual(values = c("black","orange"),
                     labels= c("Absence","Presence")) + 
  theme_bw()

# 2. Build mesh ####

mesh <- fm_mesh_2d(loc = pcod_sf, # instead of boundary, we give points
                  cutoff = 2,
                  max.edge = c(7,20), # The largest allowed triangle edge length.
                  offset = c(5,50)) # The automatic extension distance

ggplot() + 
  gg(mesh) +
  geom_sf(data= pcod_sf, aes(color = factor(present)), size = 0.5) +
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

formula <-  present ~ Intercept + 
  Eff.depth + Eff.slope + Eff.tri + 
  Eff.space

lik1 <- bru_obs(formula = formula,
                data = pcod_sf,
                family = "binomial")

m1 <- bru(cmp, 
          lik1)

summary(m1)

plot(spde.posterior(m1, "Eff.space", what = "range")) +
plot(spde.posterior(m1, "Eff.space", what = "log.variance"))

## 5.1. Predict ####

new_data <- fm_pixels(mesh, 
                      mask = T, 
                      dims = c(120, 100)) %>% 
  st_set_crs(st_crs(pcod_sf))

new_rs <- rast(st_coordinates(new_data), type = "xyz")

pred1 <- predict(m1, 
                 new_data, 
                 ~ plogis(Intercept + Eff.depth_rw + Eff.space))


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
    depth_new = seq(min(depth$depth_group[], na.rm = T), 
                        max(depth$depth_group[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.depth_rw_eval(depth_new)) 

ggplot(depth.pred) +
  geom_line(aes(depth_new, q0.5)) +
  geom_ribbon(aes(depth_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

m1$summary.random$Eff.depth_rw %>%
  ggplot() + geom_line(aes(ID,mean)) +
  geom_ribbon(aes(ID,
                  ymin = `0.025quant`,
                  ymax = `0.975quant`),
              alpha = 0.5)
