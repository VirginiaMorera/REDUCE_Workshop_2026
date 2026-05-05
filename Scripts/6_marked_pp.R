library(inlabru)
library(INLA)

data("gorillas_sf")

nests <- gorillas_sf$nests

#simulating a mark
set.seed(42)
n <- nrow(nests)
nests$mark <- 10 + rnorm(n, mean = 0, sd = 2) 

mesh <- gorillas_sf$mesh

#Model
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
  Eff.elev(covars_s$elev, model = "linear") + 
  Eff.slope(covars_s$slopeangle, model = "linear")

lik1 <- bru_obs(formula = geometry ~ -1 + inter_point + point_field + 
                  Eff.elev + Eff.slope, 
                family = "cp",
                data = nests,
                domain =  list(geometry = mesh))

lik2 <- bru_obs(formula = mark ~ -1 + inter_mark + mark_field + scale*point_field + 
                  Eff.elev + Eff.slope, 
                family = "gaussian",
                data = nests,
                domain =  list(geometry = mesh))

m1 <- bru(cmp, lik1, lik2, 
          options = list(control.inla = list(int.strategy = "eb")))

summary(m1)

newdf <- fm_pixels(mesh,
                   dims = c(100, 100),
                   mask = gorillas_sf$boundary,
                   format = "sf")

pred_points <- predict(
  m1, 
  newdf, 
  formula = ~ exp(inter_point + point_field))

pred_marks <- predict(
  m1, 
  newdf, 
  formula = ~ inter_mark + scale*point_field + mark_field)

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

