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

# 1. Obtaining data ####

data("gorillas_sf")

nests <- gorillas_sf$nests

#simulating a mark
set.seed(42)
n <- nrow(nests)
nests$mark <- 10 + rnorm(n, mean = 0, sd = 2) 

mesh <- gorillas_sf$mesh

# 2. Model ####
matern_points <- inla.spde2.pcmatern(mesh,
                                     prior.range = c(2, 0.1),
                                     prior.sigma = c(0.1, 0.1))

matern_marks <- inla.spde2.pcmatern(mesh,
                                    prior.range = c(2, 0.1),
                                    prior.sigma = c(0.1, 0.1))
cmp <- ~ -1 +
  point_field(geometry, model = matern_points) +
  mark_field(geometry, model = matern_marks) +
  inter_point(1) + inter_mark(1) + scale(1) + 
  Eff.elev_points(covars_s$elev, model = "linear") + 
  Eff.elev_marks(covars_s$elev, model = "linear") + 
  Eff.slope_points(covars_s$slopeangle, model = "linear") + 
  Eff.slope_marks(covars_s$slopeangle, model = "linear")

lik1 <- bru_obs(formula = geometry ~ -1 + inter_point + point_field + 
                  Eff.elev_points + Eff.slope_points, 
                family = "cp",
                data = nests,
                domain =  list(geometry = mesh))

lik2 <- bru_obs(formula = mark ~ -1 + inter_mark + mark_field + scale*point_field + 
                  Eff.elev_marks + Eff.slope_marks, 
                family = "gaussian",
                data = nests,
                domain =  list(geometry = mesh))

m1 <- bru(cmp, lik1, lik2, 
          options = list(control.inla = list(int.strategy = "eb")))

summary(m1)

# 3. Predict ####

newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = gorillas_sf$boundary,
                   format = "sf")

pred_points <- predict(
  m1, 
  newdf, 
  formula = ~ exp(inter_point + point_field + Eff.elev_points + Eff.slope_points))

pred_marks <- predict(
  m1, 
  newdf, 
  formula = ~ inter_mark + scale*point_field + mark_field + Eff.elev_marks + Eff.slope_marks)

ggplot() + 
  gg(data = pred_points, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, col = "red", size = 0.5) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() + 

ggplot() + 
  gg(data = pred_marks, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = nests, aes(col = mark), size = 1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(option = "B") + 
  theme_bw() 

# 4. Random spatial effect ####
rdm_points <- predict(
  m1, 
  newdf, 
  formula = ~ point_field)

rdm_marks <- predict(
  m1, 
  newdf, 
  formula = ~ mark_field)

ggplot() + 
  gg(data = rdm_points, aes(fill = q0.5), geom = "tile") +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm_points$q0.5)) * c(-1, 1)) + 
  theme_bw() + 
  
ggplot() + 
  gg(data = rdm_marks, aes(fill = q0.5), geom = "tile") +
  scale_fill_distiller(palette = 'RdBu', 
                       limit = max(abs(rdm_marks$q0.5)) * c(-1, 1)) + 
  theme_bw() 

elev.points.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elev_points_eval(elevation_new)) 

elev.marks.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    elevation_new = seq(min(covars_s$elev[], na.rm = T), 
                        max(covars_s$elev[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.elev_marks_eval(elevation_new)) 

ggplot(elev.points.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() +

ggplot(elev.marks.pred) +
  geom_line(aes(elevation_new, q0.5)) +
  geom_ribbon(aes(elevation_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

slope.points.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    slope_new = seq(min(covars_s$slopeangle[], na.rm = T), 
                        max(covars_s$slopeangle[], na.rm = T), 
                        length.out = 100)),
  formula = ~ Eff.slope_points_eval(slope_new)) 

slope.marks.pred <- predict(
  m1,
  n.samples = 100,
  newdata = data.frame(
    slope_new = seq(min(covars_s$slopeangle[], na.rm = T), 
                    max(covars_s$slopeangle[], na.rm = T), 
                    length.out = 100)),
  formula = ~ Eff.slope_marks_eval(slope_new)) 

ggplot(slope.points.pred) +
  geom_line(aes(slope_new, q0.5)) +
  geom_ribbon(aes(slope_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw() +
  
  ggplot(slope.marks.pred) +
  geom_line(aes(slope_new, q0.5)) +
  geom_ribbon(aes(slope_new,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) + 
  theme_bw()

bru_names(m1)

plot(m1, "scale")