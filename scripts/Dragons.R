# Example taken from
# https://ourcodingclub.github.io/tutorials/mixed-models/

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


load("~/Documents/GitHub/MultilevelModels/data/dragons.RData")

glimpse(dragons)

ggplot( dragons, mapping= aes( x = testScore)) +
  geom_histogram() + ggtitle("Test Score histogram") 

#ScaledDragons

sdragons= dragons %>% select(testScore, bodyLength,
                             mountainRange, site) %>% 
          mutate( across(  .cols = bodyLength, .fns = scale    ))


# Pooled data. Ignoring sited and mountain ranges.

M1= lm(  testScore ~ bodyLength, data = sdragons   )

summary(M1)

ggplot( sdragons, aes(x= bodyLength, y=testScore))+
  geom_point() +
  geom_smooth( method= "lm")

# Size is affecting test score!
#lets see the residuals 
r2(M1)
plot(check_normality(M1))
plot(check_heteroscedasticity(M1))
check_model(M1, check= "normality")
check_model(M1, check= "qq")
check_model(M1, check= "outliers")


ggplot( sdragons, aes( x= mountainRange, y=testScore, fill= mountainRange ))+
  geom_boxplot()
ggplot( sdragons, aes( x= site, y=testScore, fill= site ))+
  geom_boxplot()

#Groups behave differently
ggplot(sdragons, aes(x = bodyLength, y = testScore, colour = mountainRange)) +
  geom_point(size = 2) +
  theme_classic() 

# Plot according to group

# Considering mountainRange
ggplot( sdragons, aes(x= bodyLength, y=testScore))+
  geom_point() +
  facet_wrap( ~ mountainRange  ) +
  geom_smooth( method= "lm")

# Considering site
ggplot( sdragons, aes(x= bodyLength, y=testScore))+
  geom_point() +
  facet_wrap( ~ site  ) +
  geom_smooth( method= "lm")

# 8 mountaings, 3 sites and 2 parameters for each = 48 parameters
# Regressions are calculated independently per group. We are not considering 
# pooling here.

M2 = lm( testScore ~ bodyLength  + mountainRange  , data= sdragons )
summary(M2)
#bodylength is not significant now. M2 is testing difference in scores
# between mountain ranges. This is not our interest. Our interest is to see
# if body length affects test scores and to account or control variability that comes
# from mountain ranges and sites.

# This makes us think of random effects.

#To know if somethins is a fixed or random effect see: https://dynamicecology.wordpress.com/2015/11/04/is-it-a-fixed-or-random-effect/

#It is hard at the beginning (and maybe always) to differentiate properly
#if something should be a fixed or random effect.

M3 <-  lmer(testScore ~ bodyLength + (1| mountainRange)    , data= sdragons   )

summary(M3)

# Differences in mountain ranges account for 60% of variability after the variance
# explained by the fixed effects. Now it is clear that bodylength is not the one
# affecting intelligence!

check_model(M3)
icc(M3)

#Type of random effects

# There are “hierarchical linear models” (HLMs) or “multilevel models” out there, but while all HLMs
# are mixed models, not all mixed models are hierarchical. That’s because you can have crossed (or partially crossed)
# random factors that do not represent levels in a hierarchy.

# A factor is fully crossed if all subjects have experienced all levels of that effect.
# A factor is nested if it is inside another factor, it forms a hierarchy. 

# Data collected within sites may be correlated, hence we should include it as an additional
# random effect.

# The nesting of the site within the mountain range is implicit. They are meaningless
# w/o being assinged to a mountain range. We need a variable that denotes
# site per mountain range

sdragons <- within(sdragons, sample <- factor(mountainRange:site))

# Summing up: 
# nested random effects, the factor appears ONLY within a particular level of another factor 
# crossed effects, a given factor appears in more than one level of another factor 

# Wrong model (we have discussed why already)
mixed.WRONG <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|site), data = sdragons)

#both factors taken into account with random intercepts 
# taking nesting into account
# the notation implies independence bw mountainRange and sample site

M4 = lmer(testScore ~ bodyLength + (1|mountainRange) + (1|sample), data = sdragons)  
summary(M4)

#scaled data 
ggplot(sdragons, aes(x = bodyLength, y = testScore, colour = site)) +
    facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(sdragons, pred = predict(M4)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels

#non scaled data
ggplot(dragons, aes(x = bodyLength, y = testScore, colour = site)) +
  facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(dragons, pred = predict(M4)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels


#Including random slopes now

M5 <- lmer(testScore ~ bodyLength + (1 + bodyLength|mountainRange/site), data = sdragons) 
summary(M5)
check_model(M5)



#non scaled data
ggplot(dragons, aes(x = bodyLength, y = testScore, colour = site)) +
  facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(dragons, pred = predict(M5)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels


#body length does not affect how smart a dragon is, we can pick a small dragon to train

#---- brms

b1 <- brm(testScore ~ bodyLength + (1 + bodyLength|mountainRange/site), data = sdragons,
          iter=10000, warmup = 1000, chains=1, cores=8, seed=12) 

summary(b1)
pp_check(sdragons$testScore,  )




