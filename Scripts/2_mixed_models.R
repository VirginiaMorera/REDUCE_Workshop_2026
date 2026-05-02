# 0. Housekeeping ####
rm(list = ls())

library(INLA)
library(inlabru)
library(tidyverse)
library(ggplot2)
library(patchwork)


# 1. Very simple example

## 1.1 Generate data 

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

## 1.2 Model components and run

cmp1 = ~ Intercept(1) + 
  beta_1(x, model = "linear") +
  u(j, model = "iid")

lik1 = bru_obs(formula = y ~ Intercept + beta_1 + u,
              family = "gaussian",
              data = df)

m1 <- bru(cmp1, 
          lik1, 
          options = list(control.compute = list(config = TRUE)))

summary(m1)


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


## 1.3 predict 

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
