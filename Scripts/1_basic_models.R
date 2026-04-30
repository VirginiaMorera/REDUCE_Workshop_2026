# 0. Housekeeping ####
rm(list = ls())

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)
# load some libraries to generate nice plots
library(scico)


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

# 1.2 Run model: the simple way #### 
m0 <- bru(y ~ x, data = df)
summary(m0)

# 1.3 Run model: the "slightly" more complicated way ####

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

# 1.4 Predict ####

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

# 1.5 Priors (don't panic) ####

inla.priors.used(m0.1)

# Actually change the priors for the betas
cmp2 = ~ Intercept(1, mean.linear = 1, prec.linear = 0.1) +
  beta_1(x, model = "linear", mean.linear = 1, prec.linear = 0.1)

m0.2 = bru(cmp2, lik)

summary(m0.1)
summary(m0.2)

# 1.6 Visualise posterior marginal distributions of fixed effects
bru_names(m0.1)

plot(m0.1, "Intercept") +
plot(m0.1, "beta_1")


plot(m0.2, "Intercept") +
  plot(m0.2, "beta_1")
