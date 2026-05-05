# 0. Housekeeping ####
rm(list = ls())
load("Environments/3_basic_spatial_models.RData")

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
  geom_sf(data = gorillas_sf$plotsample$plots, fill = NA) +
  theme_bw()

plotsamples <- gorillas_sf$plotsample$counts %>% 
  mutate(PA = if_else(count == 0, 0, 1))

# 2. Covars ####
covars <- rast("Data/covars.grd")
plot(covars)

covars_s <- scale(covars)
plot(covars_s)

ggplot() + 
  geom_spatraster(data = covars_s, aes(fill = elev)) +
  scale_fill_viridis_c() +
  geom_sf(data = plotsamples, aes(col = as.factor(PA)), size = 3) +
  theme_bw()

# 3. Mesh ####

mesh <- gorillas_sf$mesh

ggplot() + 
  gg(mesh) + 
  geom_sf(data = gorillas_sf$boundary, col = "red", fill = NA) + 
  geom_sf(data = plotsamples, aes(col = as.factor(PA)), size = 3) +
  theme_bw()

# 4. Run binomial model ####
cmp <- ~  Intercept(1)  +  
  Eff.elevation(covars_s$elev, model = "linear") + 
  Eff.slope(covars_s$slopeangle, model = "linear") + 
  Eff.water(covars_s$waterdist, model = "linear")

lik1 = bru_obs(formula = PA ~ Intercept + Eff.elevation,
               family = "binomial",
               data = plotsamples, 
               domain =  list(geometry = mesh))

m1 <- bru(cmp, 
          lik1) 

summary(m1)

## 4.1 Set up dataset for predictions ####
newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = gorillas_sf$boundary,
                   format = "sf")

## 4.2 Obtain prediction ####
pred1 <- predict(m1, newdata = newdf, 
                 ~ plogis(Intercept + Eff.elevation), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = factor(PA))) +
  scale_fill_viridis_c() +
  theme_bw() 

## 4.3 Evaluate covariate effect ####
elev.pred1 <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred1) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

# 5. Run Poisson model ####

lik2 = bru_obs(formula = count ~ Intercept + Eff.elevation,
               family = "poisson",
               data = plotsamples, 
               domain =  list(geometry = mesh))

m2 <- bru(cmp, 
          lik2) 

summary(m2)

## 5.2 Obtain prediction ####
pred2 <- predict(m2, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = plotsamples, aes(col = count)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c(option = "B") + 
  theme_bw() 

## 5.3 Evaluate covariate effect ####
elev.pred2 <- predict(
  m2,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred2) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

# 6. Run log gaussian cox process model ####

## 6.1 Plot data ####

nests <- gorillas_sf$nests

ggplot() + 
  gg(mesh) +
  geom_sf(data = nests, alpha = 0.2) + 
  theme_bw()

## 6.2 Run model ####
lik3 = bru_obs(formula = geometry ~ Intercept + Eff.elevation,
               family = "cp",
               data = nests, 
               # samplers = gorillas_sf$plotsample$plots,
               domain =  list(geometry = mesh))

m3 <- bru(cmp, 
          lik3) 

summary(m3)

## 6.3 Obtain prediction ####
pred3 <- predict(m3, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred3, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red") +
  scale_fill_viridis_c() +
  theme_bw() 

## 6.4 Evaluate covariate effect ####
elev.pred3 <- predict(
  m3,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred3) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

# 7. Sampled lgcp ####

nests <- gorillas_sf$nests

sampled_nests <- st_filter(
  nests, gorillas_sf$plotsample$plots)

ggplot() + 
  geom_sf(data = nests, alpha = 0.5) + 
  geom_sf(data = gorillas_sf$plotsample$plots, fill = NA) + 
  geom_sf(data = sampled_nests, col = "red") + 
  theme_bw()

## 7.1 Run model ####
lik4 = bru_obs(formula = geometry ~ Intercept + Eff.elevation,
               family = "cp",
               data = sampled_nests, 
               samplers = gorillas_sf$plotsample$plots,
               domain =  list(geometry = mesh))

m4 <- bru(cmp, 
          lik4) 

summary(m4)

## 7.2 Obtain prediction ####
pred4 <- predict(m4, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred4, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = sampled_nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() +  

ggplot() + 
  gg(data = pred3, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() 

## 7.3 Evaluate covariate effect ####
elev.pred4 <- predict(
  m4,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

ggplot(elev.pred4) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()


m1$summary.fixed
m2$summary.fixed
m3$summary.fixed
m4$summary.fixed

ggplot(elev.pred1) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 
  
ggplot(elev.pred2) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 
  
ggplot(elev.pred3) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 
  
ggplot(elev.pred4) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 
  
plot_layout(nrow = 2)
  
# 8. Point pattern with more covariates ####
layerCor(covars_s, "cor")$correlation

lik5 <- bru_obs(formula = geometry ~ Intercept + Eff.elevation + 
                  Eff.slope + Eff.water,
                family = "cp",
                data = sampled_nests, 
                samplers = gorillas_sf$plotsample$plots,
                domain =  list(geometry = mesh))

m5 <- bru(cmp, 
          lik5)

summary(m5)

pred5 <- predict(m5, newdata = newdf, 
                 ~ exp(Intercept + Eff.elevation + 
                         Eff.slope + Eff.water), 
                 n.samples = 100)

ggplot() + 
  gg(data = pred5, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() +  
  
ggplot() + 
  gg(data = pred4, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() 

elev.pred5 <- predict(
  m5,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elevation_eval(elevation_new)) 

slope.pred5 <- predict(
  m5,
  n.samples = 100,
  newdata = data.frame(
    slope_new = seq(min(covars_s$slopeangle[], na.rm = T), 
                        max(covars_s$slopeangle[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.slope_eval(slope_new)) 

water.pred5 <- predict(
  m5,
  n.samples = 100,
  newdata = data.frame(
    water_new = seq(min(covars_s$waterdist[], na.rm = T), 
                        max(covars_s$waterdist[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.water_eval(water_new)) 

ggplot(elev.pred5) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 

ggplot(slope.pred5) +
  geom_line(aes(slope_new, q0.5)) +
  geom_ribbon(aes(slope_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() + 
  
ggplot(water.pred5) +
  geom_line(aes(water_new, q0.5)) +
  geom_ribbon(aes(water_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()


