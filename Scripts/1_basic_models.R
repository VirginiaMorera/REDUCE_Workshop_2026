# 0. Housekeeping ####
rm(list = ls())

library(INLA)
library(inlabru)
library(tidyverse)
library(ggplot2)
library(patchwork)


# 1. Very simple example ####

## 1.1 Simulate data ####

# set seed for reproducibility
set.seed(1234)

beta = c(2,0.5)
sd_error = 0.1
n = 100

x = rnorm(n)

y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)

ggplot(df) + 
  geom_point(aes(x = x, y = y)) + 
  theme_bw()

## 1.2 Run model: the simple way #### 
m0 <- bru(y ~ x, data = df)
summary(m0)

## 1.3 Run model: the "slightly" more complicated way ####

# components of the response
cmp <-  ~ Intercept(1) + beta_1(x, model = "linear")

# formula 
form <-  y ~ Intercept + beta_1

# likelihood
lik <-  bru_obs(formula = form,
              family = "gaussian",
              data = df)
# run the model
m0.1 <-  bru(cmp, lik)

summary(m0)
summary(m0.1)

## 1.4 Predict ####

new_data = data.frame(x = rnorm(100))

pred = predict(m0.1, new_data, ~ Intercept + beta_1,
               n.samples = 1000)


ggplot() +
  geom_point(data = df, aes(x = x, y = y), alpha = 0.3) +
  geom_point(data = pred, aes(x = x, y = q0.5), col = "red") +
  geom_line(data = pred, aes(x = x, y = q0.5)) +
  geom_line(data = pred, aes(x = x, y = q0.025), linetype = "dashed") +
  geom_line(data = pred, aes(x = x, y = q0.975), linetype = "dashed") +
  labs(x = "Covariate", y = "Observations") + 
  theme_bw()

## 1.5 Priors (don't panic) ####

inla.priors.used(m0.1)

# Actually change the priors for the betas
cmp2 = ~ Intercept(1, mean.linear = 1, prec.linear = 0.1) +
  beta_1(x, model = "linear", mean.linear = 1, prec.linear = 0.1)

m0.2 = bru(cmp2, lik)

summary(m0.1)
summary(m0.2)

## 1.6 Visualise posterior marginal distributions of fixed effects
bru_names(m0.1)

plot(m0.1, "Intercept") +
plot(m0.1, "beta_1")


plot(m0.2, "Intercept") +
  plot(m0.2, "beta_1")

## 1.7 Compare with normal LM ####

summary(lm(y ~ x, data = df))
confint(lm(y ~ x, data = df))

# 2. Real data example ####

## 2.1 Explore data ####
data(penguins) 
head(penguins)

penguins <- penguins %>% 
  drop_na() 

ggplot(penguins) + 
  geom_point(aes(x = flipper_len, y = body_mass)) + 
  geom_smooth(aes(x = flipper_len, y = body_mass), method = "lm") + 
  theme_bw()

## 2.2 Run model ####

cmp <- ~ Intercept(1) + flipper(flipper_len, model = "linear")

lik <- bru_obs(formula = body_mass ~ Intercept + flipper, 
               family = "Gaussian", 
               data = penguins)

m1 <- bru(cmp, lik) 
summary(m1)

## 2.3 compare with normal lm ####
summary(lm(body_mass ~ flipper_len, data = penguins))

# !!!!

## 2.4 Scaled variables ####

penguins <- penguins %>% 
  mutate(bm = scale(body_mass),
         fl = scale(flipper_len))


cmp2 <- ~ Intercept(1) + flipper(fl, model = "linear")
lik2 <- bru_obs(formula = bm ~ Intercept + flipper, 
               family = "Gaussian", 
               data = penguins)

m1.2 <- bru(cmp2, lik2) 
summary(m1.2)

## 2.5 compare with normal lm ####
summary(lm(bm ~ fl, data = penguins))
confint(lm(bm ~ fl, data = penguins))

## 2.6 Predict ####

new_data <- data.frame(
  fl = runif(100, min = min(penguins$fl), max = max(penguins$fl))
)
preds = predict(m1.2, new_data, 
                ~ Intercept + flipper, 
                samples = 1000)

ggplot() +
  geom_point(data = penguins, aes(x = fl, y = bm), alpha = 0.3) +
  geom_point(data = preds, aes(x = fl, y = q0.5), col = "red") +
  geom_line(data = preds, aes(x = fl, y = q0.5)) +
  geom_line(data = preds, aes(x = fl, y = q0.025), linetype = "dashed") +
  geom_line(data = preds, aes(x = fl, y = q0.975), linetype = "dashed") +
  labs(x = "Covariate", y = "Observations") + 
  theme_bw()

