# Example taken from
# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multilevel-models.html#multilevel-posterior-predictions

library(mvtnorm)
library(purrr)
library(tidyr)
library(stringr)
library(janitor)
library(tidybayes)
library(forcats)
library(rstan) 
library(brms)
library(ggplot2)
library(dplyr)
library(loo)
library(nlme)
library(lme4)
library(GGally)
library(corrplot)
library(tidyverse)
library(bayesplot)
library(wesanderson)
library(cmdstanr)
library(performance)
library(see)
library(ggrepel)
library(qqplotr)


setwd("~/Documents/GitHub/MultilevelModels/scripts")
PATH_TO_CMDSTAN= "~/cmdstan"
set_cmdstan_path(PATH_TO_CMDSTAN)
cmdstan_path()
cmdstan_version()
options(mc.cores = parallel::detectCores())

file <- "~/Documents/GitHub/MultilevelModels/data/reedfrogs.csv"
frogs=read.table( file,sep=";" , header= TRUE)

frogs %>% glimpse()

frogs %>% 
  ggplot( aes( x=size, y=surv, fill=surv ) )+
  geom_boxplot()

frogs %>% 
  ggplot( aes( x=pred, y=surv, fill=pred ) )+
  geom_boxplot()



#Make a cluster

frogs = frogs %>% 
        mutate( tank= 1:nrow(frogs)   )

# Unpooled model. Each tank has an intercept.
# surv ~ bin(ni, pi )
# logit(pi)= alpha_tank
# alpha_tank ~Normal(0,5)



b1 = brm( data = frogs, family = binomial, 
          surv | trials(density) ~ 0+factor(tank),
          prior(normal(0, 5), class = b),
          iter = 4000, warmup = 500, chains = 4, cores = 8,
          seed = 12
          )


# Pooled model. Each tank has an intercept.
# surv ~ bin(ni, pi )
# logit(pi)= alpha_tank
# alpha_tank ~Normal( alpha, sigma)
# alpha ~N(0,1)
# sigma ~C+(0,3)

b2 = brm( data = frogs, family = binomial, 
          surv | trials(density) ~ 1+(1|tank),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 3), class = sd)),
          iter = 4000, warmup = 500, chains = 4, cores = 8,
          seed = 12
        )


# WAIC comparisons
b1 <- add_criterion(b1, "waic")
b2 <- add_criterion(b2, "waic")

w=loo_compare(b1, b2, criterion= "waic"   )

print(w, simplify=F)

#posterior samples of b2
postb2 = posterior_samples(b2, add_chain=T)

ppred=coef(b2)$tank[,,] %>% 
  as_tibble() %>% 
  bind_cols(frogs) %>% 
  mutate(  ppred= inv_logit_scaled(Estimate)    )


#this is a nice plot
ppred %>% 
ggplot(aes(x = tank)) +
  geom_hline(yintercept = inv_logit_scaled(median(postb2$b_Intercept)), linetype = 2, size = 1/4)+
  geom_vline(xintercept = c(16.5, 32.5), size = 1/4) +
  geom_point(aes(y = propsurv), color = "orange2") +
  geom_point(aes(y = ppred), shape = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(1, 16, 32, 48)) +
  labs(title    = "Multilevel shrinkage!",
       subtitle = "The empirical proportions are in orange while the model-\nimplied proportions are the black circles. The dashed line is\nthe model-implied average survival proportion.") +
  annotate("text", x = c(8, 16 + 8, 32 + 8), y = 0, 
           label = c("small tanks", "medium tanks", "large tanks")) 


ggplot(data = postb2, 
       aes(x = rnorm(n    = nrow(postb2), 
                     mean = b_Intercept, 
                     sd   = sd_tank__Intercept) %>% 
             inv_logit_scaled())) +
  geom_density(size = 1, fill = "blue") +
  ggtitle("Probability of survival") 
  

#summary of b2
print(b2)


# 3 perspectives
# Complete Pooling (single alpha model)
# No pooling (single level alpha_tank model)
# Partial Pooling (multilevel model alpha_tank ~N(alpha, sigma))


# Idea: we will simulate some tadpole data. We will know the true per pond
# survival probabilities. We can compare no pooling estimates to partial pooling
# estimates 

#Simulation

# surv ~ bin(ni, pi )
# logit(pi)= alpha_tank
# alpha_tank ~Normal( alpha, sigma)
# alpha ~N(0,1)
# sigma ~C+(0,3)



a=1.4
sigma=1.5
nponds=60

set.seed(1)

fsim= tibble(
      pond= 1:nponds,
      ni=rep(c(5,10,25,35), each= nponds/4) %>%  as.integer(),
      true_a = rnorm( nponds, mean=a, sd=sigma )
      )


#simulate survivors

fsim= fsim %>% 
  mutate( si=rbinom(n=n(), prob=inv_logit_scaled(true_a),size=ni))

#no pooling estimates

fsim = fsim %>% 
      mutate( p_nopool= si/ni)


#partial pooling estimates

b3 = brm( data= fsim, family=binomial,
          si | trials(ni) ~1+(1|pond),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 3), class = sd)),
          iter=10000, warmup=1000, chains=1, cores=8, seed=12
      )

print(b3)

p_partpool <- 
  coef(b3)$pond[, , ] %>% 
  as_tibble() %>%
  transmute(p_partpool = inv_logit_scaled(Estimate))

fsim <- 
  fsim %>%
  bind_cols(p_partpool) %>% 
  mutate(p_true         = inv_logit_scaled(true_a)) %>%
  mutate(nopool_error   = abs(p_nopool   - p_true),
         partpool_error = abs(p_partpool - p_true))

fsim %>% 
  glimpse()

#compare no pooling vs pooling
fsim %>%
  select(ni, nopool_error:partpool_error) %>%
  gather(key, value, -ni) %>%
  group_by(key) %>%
  summarise(mean_error   = mean(value) %>% round(digits = 3),
            median_error = median(value) %>% round(digits = 3))




  



