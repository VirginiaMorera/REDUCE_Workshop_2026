# 0. Housekeeping ####
rm(list = ls())

library(INLA)
library(inlabru)
library(tidyverse)
library(ggplot2)
library(patchwork)

bru_options_set(control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       cpo = TRUE))

# 1. Linear mixed model ####

## 1.1 Generate data ####

set.seed(12)
beta <-  c(1.5,1)
sd_error <-  1
tau_group <-  1
n <-  100
n.groups  <-  5
x <-  rnorm(n)
v <-  rnorm(n.groups, sd = tau_group^{-1/2})
y <-  beta[1] + beta[2] * x + rnorm(n, sd = sd_error) +
  rep(v, each = 20)
df <-  data.frame(y = y, x = x, j = rep(1:5, each = 20))

ggplot(df) +
  geom_point(aes(x = x, colour = factor(j), y = y)) +
  theme_classic() +
  scale_colour_discrete("Group")

## 1.2 Model components and run ####

cmp1 = ~ Intercept(1) + 
  beta_1(x, model = "linear") +
  u(j, model = "iid")

lik1 = bru_obs(formula = y ~ Intercept + beta_1 + u,
              family = "gaussian",
              data = df)

m1 <- bru(cmp1, 
          lik1)

summary(m1)

## 1.3 Changing priors ####

inla.priors.used(m1)

cmp1.2 = ~ Intercept(1) + 
  beta_1(x, model = "linear") +
  u(j, model = "iid", 
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.01)
      )
    ))

m1.2 <- bru(cmp1.2, 
            lik1)

summary(m1)
summary(m1.2)

## 1.4 Evaluate random effect ####

samples <- inla.posterior.sample(n = 100, result = m1)

hyper_df <- map_dfr(samples, ~ {
  tibble(
    name  = names(.x$hyperpar),
    value = as.numeric(.x$hyperpar)
  )
}, .id = "sample_id")

hyper_df <- hyper_df %>%
  pivot_wider(names_from = name, values_from = value)

hyper_df <- hyper_df %>% 
  mutate(var_data = 1 / `Precision for the Gaussian observations`, 
         var_rdm = 1/ `Precision for u`, 
         icc = var_rdm / (var_rdm + var_data))         

icc_sum <-  hyper_df %>%
  summarise(mean = mean(icc), 
            sd = sd(icc), 
            q025 = quantile(icc, 0.025),
            q050 = quantile(icc, 0.5),
            q975 = quantile(icc, 0.975))

## 1.5 predict #### 

new_data <- data.frame(
  x = rnorm(100), 
  j = rep(1:5, each = 20)
)

pred  <-  predict(m1, 
                  new_data, 
                  formula = ~ Intercept + beta_1 + u, 
                  n.samples = 1000)

ggplot(data = pred, aes(x = x,y = mean, color = factor(j))) +
  geom_point() + 
  geom_line()+
  geom_ribbon(aes(x, ymin = q0.025, ymax = q0.975, fill = factor(j)), alpha = 0.2) +
  geom_point(data = df, aes(x = x,y = y), col = "gray", alpha = 0.5) +
  # facet_wrap(~j) + 
  theme_bw()

# 2. Generalised linear model ####

## 2.1 Generate data ####

set.seed(123)
n <-  100
beta <- c(1,1)
x <-  rnorm(n)
lambda <-  exp(beta[1] + beta[2] * x)
y <-  rpois(n, lambda = lambda)
df <-  data.frame(y = y, x = x)


ggplot(df) + 
  geom_point(aes(x = x, y = y), alpha = 0.5) + 
  geom_smooth(aes(x = x, y = y), 
              method = "glm", 
              method.args = list(family = "poisson"))

## 2.2 Model components and run ####

cmp2 <- ~ Intercept(1) + beta_1(x, model = "linear")

lik2 <- bru_obs(formula = y ~ Intercept + beta_1,
              family = "poisson",
              data = df)

m2 <- bru(cmp2, lik2)

summary(m2)

## 2.3 Changing priors ####

inla.priors.used(m2)

m2.2 <- bru(cmp2, lik2,  
            options = list(
              control.fixed = list(
                mean = 0,
                prec = 100)
            ))

summary(m2)
summary(m2.2)

## 2.4 Predict ####

new_data <- data.frame(
  x = runif(100, min(df$x), max(df$x))
)

pred <- predict(m2, new_data, 
                ~ exp(Intercept + beta_1), 
                nsamples = 1000)

ggplot() +
  geom_point(data = df, aes(x = x, y = y), alpha = 0.3) +
  geom_point(data = pred, aes(x = x, y = q0.5), col = "red") +
  geom_line(data = pred, aes(x = x, y = q0.5)) +
  geom_line(data = pred, aes(x = x, y = q0.025), linetype = "dashed") +
  geom_line(data = pred, aes(x = x, y = q0.975), linetype = "dashed") +
  labs(x = "Covariate", y = "Observations") + 
  theme_bw()


# 3. GAMs the inlabru way ####

## 3.1 Generate data ####

set.seed(123)
n <-  100
x <-  rnorm(n)
eta <-  (1 + cos(x))
y <-  rnorm(n, mean = eta, sd = 0.5)
df <-  data.frame(y = y,
                  x = inla.group(x)) # equidistant x's

ggplot(df) + 
  geom_point(aes(x = x, y = y), alpha = 0.3) + 
  geom_smooth(aes(x = x, y = y)) + 
  theme_bw()

## 3.2 Model components and run ####

cmp3 <- ~ Intercept(1) +
  smooth(x, model = "rw1")

lik3 <- bru_obs(formula = y ~.,
              family = "gaussian",
              data = df)

m3 <- bru(cmp3, lik3)

summary(m3)

smooth <- m3$summary.random$smooth

ggplot(smooth) +
  geom_ribbon(aes(x = ID, 
                  ymin = `0.025quant`, 
                  ymax = `0.975quant`), 
              alpha = 0.5) +
  geom_line(aes(x = ID, y = mean)) +
  labs(x = "covariate", y = "") + 
  theme_bw()

## 3.3 Changing priors ####

inla.priors.used(m3)

cmp3.2 <- ~ Intercept(1) +
  smooth(x, model = "rw2", 
         hyper = list(
           prec = list(
             prior = "pc.prec",
             param = c(0.05, 0.2)
           )
         )
  )

m3.2 <- bru(cmp3.2, lik3)

smooth2 <- m3.2$summary.random$smooth

ggplot(smooth) +
  geom_ribbon(aes(x = ID, 
                  ymin = `0.025quant`, 
                  ymax = `0.975quant`), 
              alpha = 0.5) +
  geom_line(aes(x = ID, y = mean)) +
  labs(x = "covariate", y = "") + 
  theme_bw() +

ggplot(smooth2) +
  geom_ribbon(aes(x = ID, 
                  ymin = `0.025quant`, 
                  ymax = `0.975quant`), 
              alpha = 0.5) +
  geom_line(aes(x = ID, y = mean)) +
  labs(x = "covariate", y = "") + 
  theme_bw()

## 3.4 Predict ####

new_data <- data.frame(
  x = runif(100, min = min(df$x), max = max(df$y))
)

pred  <-  predict(m3, 
                  new_data, 
                  ~ Intercept + smooth, 
                  n.samples = 1000)

ggplot(pred) +
  geom_point(data = df, aes(x,y), alpha = 0.3) +
  geom_point(aes(x = x, y = q0.5), col = "red") +
  geom_ribbon(aes(x, ymin = q0.025, ymax = q0.975), 
                  alpha = 0.2) +
  geom_line(aes(x, q0.5)) +
  labs(x = "Covariate", y = "Observations") + 
  theme_bw()
